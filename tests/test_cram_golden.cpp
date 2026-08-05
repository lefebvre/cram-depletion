// Regression lock for the CRAM solver on a matrix that actually exercises
// LU fill-in.
//
// Every other solver test builds a linear decay chain: strictly triangular, so
// Eigen::SparseLU does no pivoting and produces no fill-in, and N never exceeds
// 6. Those tests are backed by analytic Bateman solutions, which is what makes
// them trustworthy -- but the closed form is an alternating sum that loses
// precision badly for long chains and requires distinct decay constants, so it
// does not scale up. That leaves the assembly and factorization path validated
// only in its easiest case.
//
// This file closes that gap from the other side: a 256-nuclide cyclic matrix
// (fill-in ratio ~2.2x) whose CRAM48 output is frozen in cram_golden_data.hpp.
//
// What this is and is not:
//   * It IS a regression lock. It will catch any change that alters solver
//     behaviour -- a mis-mapped diagonal index, a stale reused buffer, a race.
//   * It is NOT a correctness proof. It freezes what the solver does today at
//     N=256, where we have no independent reference. It cannot tell you the
//     frozen values are right. Certification has to come from a comparison
//     against an external code on real inventories -- the VERA pin replay
//     against OpenMC on the feature/openmc-vera-benchmark branch:
//     https://github.com/lefebvre/cram-depletion/tree/feature/openmc-vera-benchmark
//     That is planned but not yet landed, so until it is, treat these vectors
//     as change detection and nothing more.
//
// The two solver paths are each checked against the golden data INDEPENDENTLY,
// never against each other. cramSolve() and CramSolver share the per-pole
// assembly logic, so a bug touching both would leave them in perfect agreement
// while both drift -- a differential check cannot see that.
//
// To regenerate after an intentional behaviour change: see the PR that added
// this file. Regenerating to make a red test go green defeats the purpose --
// establish first that the change was intended.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "cram/cram.hpp"
#include "cram/cram_solver.hpp"
#include "cram_golden_data.hpp"
#include "cyclic_chain.hpp"

using cram::CramOrder;
using cram::cramSolve;
using cram::CramSolver;
using cram_test::buildCyclicChain;

namespace {

// Entries above this fraction of the largest are compared by relative error;
// everything below is compared absolutely. A flat relative tolerance is
// meaningless on trace entries -- the smallest values here sit ~50 decades below
// the largest, where no significant digits survive.
constexpr double kSignificantFrac = 1e-6;

// Measured same-platform agreement between the two solver paths is exact (0 ulp),
// and a correct refactor of the pole loop should land within a few ulp. The gap
// to a genuine defect is enormous -- a mis-mapped diagonal or stale buffer shifts
// values by 1e-3 or more, or produces NaN. So these are set loose enough to
// absorb cross-platform FP variation (FMA contraction, vectorization) that cannot
// be measured from a single machine, while still discriminating by many orders of
// magnitude. If CI proves stable across g++/clang/MSVC they can be tightened.
constexpr double kRelTol = 1e-9;
constexpr double kAbsFrac = 1e-12;

void expectMatchesGolden(const Eigen::VectorXd& got,
                         const std::array<double, cram_test::kGoldenN>& ref) {
  ASSERT_EQ(got.size(), cram_test::kGoldenN);
  double maxRef = 0.0;
  for (double r : ref) {
    maxRef = std::max(maxRef, std::abs(r));
  }
  ASSERT_GT(maxRef, 0.0);

  for (int i = 0; i < cram_test::kGoldenN; ++i) {
    const double r = ref[static_cast<std::size_t>(i)];
    if (std::abs(r) > kSignificantFrac * maxRef) {
      EXPECT_NEAR(got(i), r, kRelTol * std::abs(r)) << "index " << i;
    } else {
      EXPECT_NEAR(got(i), r, kAbsFrac * maxRef) << "index " << i << " (trace)";
    }
  }
}

}  // namespace

// The lock is only worth anything if the matrix keeps the properties it was
// built for. If someone edits cyclic_chain.hpp and flattens it back into a
// triangular chain, every other test here would still pass while silently
// testing nothing -- this one fails instead and says why.
TEST(CramGolden, MatrixHasTheStructureTheLockDependsOn) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  ASSERT_EQ(c.A.rows(), cram_test::kGoldenN);

  int above = 0;
  int below = 0;
  int missingDiagonal = cram_test::kGoldenN;
  double maxColSum = 0.0;
  for (int k = 0; k < c.A.outerSize(); ++k) {
    double colSum = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(c.A, k); it; ++it) {
      colSum += it.value();
      if (it.row() < it.col()) {
        ++above;
      } else if (it.row() > it.col()) {
        ++below;
      } else {
        --missingDiagonal;
      }
    }
    maxColSum = std::max(maxColSum, std::abs(colSum));
  }

  // Entries on both sides of the diagonal: the matrix has cycles, so the
  // factorization cannot degenerate into forward substitution.
  EXPECT_GT(above, 100) << "matrix is no longer cyclic; fill-in is not being exercised";
  EXPECT_GT(below, 100);
  // Stable sinks contribute no triplets, so (A*dt - theta*I) has to create
  // their diagonal entry -- the case a diagonal-index optimization can break.
  EXPECT_GT(missingDiagonal, 0) << "no structurally-absent diagonal left to exercise";
  // Zero column sums are what make the mass-conservation test exact.
  EXPECT_LT(maxColSum, 1e-15) << "columns no longer sum to zero; mass test is invalid";
}

TEST(CramGolden, FreeFunctionMatchesGoldenMild) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  expectMatchesGolden(cramSolve(c.A, c.n0, cram_test::kGoldenDtMild, CramOrder::CRAM48),
                      cram_test::kGoldenMild);
}

TEST(CramGolden, FreeFunctionMatchesGoldenStiff) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  expectMatchesGolden(cramSolve(c.A, c.n0, cram_test::kGoldenDtStiff, CramOrder::CRAM48),
                      cram_test::kGoldenStiff);
}

TEST(CramGolden, CachedSolverMatchesGoldenMild) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(c.A, cram_test::kGoldenDtMild);
  expectMatchesGolden(solver.apply(c.n0), cram_test::kGoldenMild);
}

TEST(CramGolden, CachedSolverMatchesGoldenStiff) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  CramSolver solver(CramOrder::CRAM48);
  solver.prepare(c.A, cram_test::kGoldenDtStiff);
  expectMatchesGolden(solver.apply(c.n0), cram_test::kGoldenStiff);
}

// --- Reference-free invariants -------------------------------------------
// These hold for any N and any dt and do not depend on the frozen vector, so
// they keep working even if the golden data is ever regenerated.

TEST(CramGolden, ConservesMass) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  const double initial = c.n0.sum();
  ASSERT_GT(initial, 0.0);
  // Every column of A sums to zero and no daughter leaves the tracked set, so
  // exp(A*dt) is exactly mass-preserving; only rounding separates the two sums.
  //
  // The bound is set from measurement, not from how tight it can be made. Worst
  // observed at the stiff dt is 1.3e-14 on MSVC and 5.1e-15 on g++. That figure
  // moves by ~10x with the generator seed and ~3x across platforms, because it
  // is rounding accumulated through 24 complex sparse solves rather than
  // anything structural -- the matrix's column sums are ~1e-17 by construction,
  // which the structural test above asserts. 1e-11 leaves roughly 750x headroom
  // while still catching a genuine conservation defect, which shows up at 1e-6
  // or worse.
  //
  // Deliberately not tuned by reseeding: the seed-to-seed spread here runs from
  // 0 to 5e-15, so picking a stream that flatters this bound would be fitting
  // the fixture to the assertion instead of the other way round.
  for (double dt : {1.0e3, cram_test::kGoldenDtMild, 1.0e7, cram_test::kGoldenDtStiff}) {
    const Eigen::VectorXd n = cramSolve(c.A, c.n0, dt, CramOrder::CRAM48);
    EXPECT_NEAR(n.sum(), initial, 1e-11 * initial) << "dt=" << dt;
  }
}

TEST(CramGolden, NoLargeNegativeInventories) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  // CRAM is a rational approximation and admits tiny negative excursions on
  // species that should be ~0, so this bounds their size rather than requiring
  // non-negativity. A sign or structure error produces negatives many orders
  // larger than this floor.
  for (double dt : {cram_test::kGoldenDtMild, cram_test::kGoldenDtStiff}) {
    const Eigen::VectorXd n = cramSolve(c.A, c.n0, dt, CramOrder::CRAM48);
    const double floorValue = -1e-15 * n.cwiseAbs().maxCoeff();
    EXPECT_GT(n.minCoeff(), floorValue) << "dt=" << dt;
  }
}

TEST(CramGolden, SemigroupProperty) {
  const auto c = buildCyclicChain(cram_test::kGoldenN, cram_test::kGoldenSeed);
  // exp(A*dt) exp(A*dt) == exp(A*2dt). CRAM only approximates the exponential,
  // so the two differ by truncation error rather than rounding alone; measured
  // worst case is ~3e-15 of the largest entry.
  for (double dt : {cram_test::kGoldenDtMild, cram_test::kGoldenDtStiff}) {
    const Eigen::VectorXd once = cramSolve(c.A, c.n0, 2.0 * dt, CramOrder::CRAM48);
    const Eigen::VectorXd first = cramSolve(c.A, c.n0, dt, CramOrder::CRAM48);
    const Eigen::VectorXd twice = cramSolve(c.A, first, dt, CramOrder::CRAM48);
    const double tol = 1e-11 * once.cwiseAbs().maxCoeff();
    for (int i = 0; i < once.size(); ++i) {
      EXPECT_NEAR(twice(i), once(i), tol) << "dt=" << dt << ", index " << i;
    }
  }
}
