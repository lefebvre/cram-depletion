#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <thread>
#include <vector>

#include "bateman.hpp"
#include "cram/cram_solver.hpp"

using namespace cram;
using cram_test::batemanLinearChain;

namespace {

Eigen::SparseMatrix<double> linearChain(const std::vector<double>& lambda) {
  const int n = static_cast<int>(lambda.size());
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    if (lambda[i] != 0.0) {
      t.emplace_back(i, i, -lambda[i]);
      if (i + 1 < n)
        t.emplace_back(i + 1, i, lambda[i]);
    }
  }
  Eigen::SparseMatrix<double> A(n, n);
  A.setFromTriplets(t.begin(), t.end());
  return A;
}

Eigen::VectorXd unitFirst(int n) {
  Eigen::VectorXd v(n);
  v.setZero();
  v(0) = 1.0;
  return v;
}

}  // namespace

TEST(CramSolver, MatchesFreeFunction) {
  auto A = linearChain({0.07, 0.013, 0.004, 0.0});
  Eigen::VectorXd n0 = unitFirst(4);
  const double dt = 50.0;

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(A, dt);
  Eigen::VectorXd cached = solver.apply(n0);
  Eigen::VectorXd direct = cramSolve(A, n0, dt, CramOrder::CRAM48);

  ASSERT_EQ(cached.size(), direct.size());
  for (int i = 0; i < cached.size(); ++i) {
    EXPECT_NEAR(cached(i), direct(i), 1e-12) << "index " << i;
  }
}

TEST(CramSolver, RepeatedApplyMarchesCorrectly) {
  // March in two equal steps and compare against a single full-dt solve and
  // against the Bateman analytic solution.
  std::vector<double> lam{0.05, 0.02, 0.0};
  auto A = linearChain(lam);
  const double dt = 30.0;

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(A, dt);

  Eigen::VectorXd y = unitFirst(3);
  y = solver.apply(y);
  y = solver.apply(y);  // two dt sub-steps

  auto ref = batemanLinearChain(lam, 2.0 * dt);
  for (std::size_t i = 0; i < ref.size(); ++i) {
    EXPECT_NEAR(y(i), ref[i], 1e-9 * std::max(1.0, std::abs(ref[i])));
  }
}

TEST(CramSolver, RePrepareSwitchesMatrix) {
  auto A1 = linearChain({0.1, 0.0});
  auto A2 = linearChain({0.5, 0.0});
  Eigen::VectorXd n0 = unitFirst(2);

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(A1, 10.0);
  Eigen::VectorXd y1 = solver.apply(n0);
  EXPECT_NEAR(y1(0), std::exp(-1.0), 1e-10);

  solver.prepare(A2, 10.0);
  Eigen::VectorXd y2 = solver.apply(n0);
  EXPECT_NEAR(y2(0), std::exp(-5.0), 1e-10);
}

TEST(CramSolver, ApplyBeforePrepareThrows) {
  CramSolver solver;
  Eigen::VectorXd n0 = unitFirst(2);
  EXPECT_THROW(solver.apply(n0), std::logic_error);
  EXPECT_FALSE(solver.prepared());
}

TEST(CramSolver, WrongSizeThrows) {
  auto A = linearChain({0.05, 0.0});
  CramSolver solver;
  solver.prepare(A, 1.0);
  EXPECT_TRUE(solver.prepared());
  EXPECT_EQ(solver.size(), 2);

  Eigen::VectorXd wrong(3);
  wrong << 1, 2, 3;
  EXPECT_THROW(solver.apply(wrong), std::invalid_argument);
}

TEST(CramSolver, NonSquareThrows) {
  Eigen::SparseMatrix<double> A(2, 3);
  CramSolver solver;
  EXPECT_THROW(solver.prepare(A, 1.0), std::invalid_argument);
  EXPECT_FALSE(solver.prepared());
}

TEST(CramSolver, OrderIsRetained) {
  EXPECT_EQ(CramSolver(CramOrder::CRAM16).order(), CramOrder::CRAM16);
  EXPECT_EQ(CramSolver(CramOrder::CRAM48).order(), CramOrder::CRAM48);
}

// --- symbolic analysis reuse ----------------------------------------------
// The many-region case: every depletion region shares the chain's topology and
// differs only in reaction rates, so re-preparing the same solver region after
// region should analyze the sparsity once and refactorize numerically each time.

TEST(CramSolver, ReusesSymbolicAnalysisWhenPatternIsUnchanged) {
  auto A1 = linearChain({0.07, 0.013, 0.004, 0.0});
  auto A2 = linearChain({0.05, 0.020, 0.009, 0.0});  // same pattern, new values

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(A1, 10.0);
  EXPECT_FALSE(solver.reusedSymbolicAnalysis()) << "nothing to reuse on the first prepare";
  solver.prepare(A2, 10.0);
  EXPECT_TRUE(solver.reusedSymbolicAnalysis());
}

TEST(CramSolver, SymbolicReuseDoesNotChangeResults) {
  auto A1 = linearChain({0.07, 0.013, 0.004, 0.0});
  auto A2 = linearChain({0.05, 0.020, 0.009, 0.0});
  Eigen::VectorXd n0 = unitFirst(4);

  // Same final matrix reached two ways: once through a solver that had already
  // analyzed a matching pattern, once through a solver seeing it cold.
  CramSolver reused(CramOrder::CRAM48);
  reused.prepare(A1, 10.0);
  reused.prepare(A2, 10.0);
  ASSERT_TRUE(reused.reusedSymbolicAnalysis());

  CramSolver fresh(CramOrder::CRAM48);
  fresh.prepare(A2, 10.0);
  ASSERT_FALSE(fresh.reusedSymbolicAnalysis());

  const Eigen::VectorXd viaReuse = reused.apply(n0);
  const Eigen::VectorXd viaFresh = fresh.apply(n0);
  ASSERT_EQ(viaReuse.size(), viaFresh.size());
  for (int i = 0; i < viaReuse.size(); ++i) {
    EXPECT_DOUBLE_EQ(viaReuse(i), viaFresh(i)) << "index " << i;
  }
}

TEST(CramSolver, FallsBackToFullAnalysisWhenSizeChanges) {
  auto small = linearChain({0.07, 0.013, 0.0});
  auto large = linearChain({0.07, 0.013, 0.004, 0.0});

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(small, 10.0);
  solver.prepare(large, 10.0);
  EXPECT_FALSE(solver.reusedSymbolicAnalysis());
  EXPECT_EQ(solver.size(), 4);
}

TEST(CramSolver, FallsBackToFullAnalysisWhenSparsityChanges) {
  // Same size, same nonzero count, but one off-diagonal moved: the pattern
  // check has to catch this, not just compare dimensions.
  std::vector<Eigen::Triplet<double>> t1{{0, 0, -0.1}, {1, 0, 0.1}, {1, 1, -0.2}, {2, 1, 0.2}};
  std::vector<Eigen::Triplet<double>> t2{{0, 0, -0.1}, {2, 0, 0.1}, {1, 1, -0.2}, {2, 1, 0.2}};
  Eigen::SparseMatrix<double> A1(3, 3);
  Eigen::SparseMatrix<double> A2(3, 3);
  A1.setFromTriplets(t1.begin(), t1.end());
  A2.setFromTriplets(t2.begin(), t2.end());

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(A1, 10.0);
  solver.prepare(A2, 10.0);
  EXPECT_FALSE(solver.reusedSymbolicAnalysis());
}

TEST(CramSolver, RepeatedRegionSweepStaysOnTheFastPath) {
  // Ten "regions": one topology, ten different rate sets. Only the first
  // should pay for the symbolic analysis.
  CramSolver solver(CramOrder::CRAM48);
  int analyses = 0;
  for (int r = 0; r < 10; ++r) {
    const double scale = 1.0 + 0.1 * r;
    auto A = linearChain({0.07 * scale, 0.013 * scale, 0.004 * scale, 0.0});
    solver.prepare(A, 10.0);
    if (!solver.reusedSymbolicAnalysis())
      ++analyses;
    EXPECT_TRUE(solver.apply(unitFirst(4)).allFinite()) << "region " << r;
  }
  EXPECT_EQ(analyses, 1);
}

// --- decomposed prepare ----------------------------------------------------

TEST(CramSolver, DecomposedPrepareMatchesPrepare) {
  auto A = linearChain({0.07, 0.013, 0.004, 0.020, 0.0});
  Eigen::VectorXd n0 = unitFirst(5);

  CramSolver whole(CramOrder::CRAM48);
  whole.prepare(A, 25.0);

  CramSolver split(CramOrder::CRAM48);
  split.beginPrepare(A, 25.0);
  ASSERT_EQ(split.poleCount(), whole.poleCount());
  for (std::size_t l = 0; l < split.poleCount(); ++l)
    split.preparePole(l);
  split.endPrepare();
  ASSERT_TRUE(split.prepared());

  const Eigen::VectorXd a = whole.apply(n0);
  const Eigen::VectorXd b = split.apply(n0);
  for (int i = 0; i < a.size(); ++i) {
    EXPECT_DOUBLE_EQ(a(i), b(i)) << "index " << i;
  }
}

TEST(CramSolver, PreparePoleOutsideBeginEndThrows) {
  CramSolver solver(CramOrder::CRAM48);
  EXPECT_THROW(solver.preparePole(0), std::logic_error);
  EXPECT_THROW(solver.endPrepare(), std::logic_error);
}

TEST(CramSolver, PreparePoleRejectsOutOfRangeIndex) {
  auto A = linearChain({0.1, 0.0});
  CramSolver solver(CramOrder::CRAM48);
  solver.beginPrepare(A, 1.0);
  EXPECT_THROW(solver.preparePole(solver.poleCount()), std::out_of_range);
}

// Drives the poles from several threads at once, which is the whole point of
// the decomposition. Under the ThreadSanitizer CI job this is what proves the
// concurrency claim in the header; without TSan it still checks that a
// thread-driven prepare produces exactly the sequential answer.
TEST(CramSolver, ConcurrentPreparePoleMatchesSequential) {
  auto A = linearChain({0.07, 0.013, 0.004, 0.020, 0.008, 0.0});
  Eigen::VectorXd n0 = unitFirst(6);
  const double dt = 3600.0;

  CramSolver sequential(CramOrder::CRAM48);
  sequential.prepare(A, dt);

  CramSolver threaded(CramOrder::CRAM48);
  threaded.beginPrepare(A, dt);

  const std::size_t poles = threaded.poleCount();
  const unsigned hw = std::max(2u, std::min(8u, std::thread::hardware_concurrency()));
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> workers;
  workers.reserve(hw);
  for (unsigned w = 0; w < hw; ++w) {
    workers.emplace_back([&] {
      for (std::size_t l = next++; l < poles; l = next++)
        threaded.preparePole(l);
    });
  }
  for (auto& t : workers)
    t.join();
  threaded.endPrepare();

  const Eigen::VectorXd a = sequential.apply(n0);
  const Eigen::VectorXd b = threaded.apply(n0);
  for (int i = 0; i < a.size(); ++i) {
    EXPECT_DOUBLE_EQ(a(i), b(i)) << "index " << i;
  }
}

// Step-by-step parity with cramSolve over a long march. Both paths perform
// the same arithmetic in the same order, so the inventories should stay
// bit-identical (or within rounding noise) across every step. If this test
// drifts, something has diverged in the per-pole factor/solve sequence.
TEST(CramSolver, LongMarchMatchesFreeFunction) {
  auto A = linearChain({0.07, 0.013, 0.004, 0.020, 0.008, 0.0});
  Eigen::VectorXd n0 = unitFirst(6);
  const double dt = 86400.0;
  const int steps = 500;

  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(A, dt);

  Eigen::VectorXd cached = n0;
  Eigen::VectorXd direct = n0;
  for (int s = 0; s < steps; ++s) {
    cached = solver.apply(cached);
    direct = cramSolve(A, direct, dt, CramOrder::CRAM48);
    ASSERT_EQ(cached.size(), direct.size());
    for (int i = 0; i < cached.size(); ++i) {
      // Use a per-step bound so we localise any drift to its first step.
      const double tol = 1e-12 * std::max(1.0, std::abs(direct(i)));
      ASSERT_NEAR(cached(i), direct(i), tol) << "step " << s << ", index " << i;
    }
  }
}
