#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
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
