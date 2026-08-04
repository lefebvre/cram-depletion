#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "bateman.hpp"
#include "cram/cram.hpp"

using namespace cram;
using cram_test::batemanLinearChain;
using cram_test::batemanTwoBranch;

namespace {

// Build the burnup matrix for a linear chain 0 -> 1 -> ... with the given decay
// constants. A terminal lambda of 0 marks a stable nuclide.
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

// Compare significant entries by relative error and near-zero entries by
// absolute error (CRAM resolves the dominant species, not values many orders
// of magnitude below them).
void expectClose(const Eigen::VectorXd& got, const std::vector<double>& ref, double relTol,
                 double absTol) {
  ASSERT_EQ(static_cast<std::size_t>(got.size()), ref.size());
  double scale = 0.0;
  for (double r : ref)
    scale = std::max(scale, std::abs(r));
  for (std::size_t i = 0; i < ref.size(); ++i) {
    if (std::abs(ref[i]) > 1e-6 * scale)
      EXPECT_NEAR(got(i), ref[i], relTol * std::abs(ref[i])) << "index " << i;
    else
      EXPECT_NEAR(got(i), ref[i], absTol) << "index " << i << " (trace)";
  }
}

constexpr double kRel48 = 1e-9;
constexpr double kRel16 = 1e-4;
constexpr double kAbs = 1e-10;

}  // namespace

// ---------------------------------------------------------------------------
TEST(Cram, ScalarExponential) {
  for (CramOrder order : {CramOrder::CRAM16, CramOrder::CRAM48}) {
    Eigen::SparseMatrix<double> A(1, 1);
    A.insert(0, 0) = -0.3;
    A.makeCompressed();
    Eigen::VectorXd n0(1);
    n0 << 5.0;
    double t = 4.0;
    Eigen::VectorXd n = cramSolve(A, n0, t, order);
    EXPECT_NEAR(n(0), 5.0 * std::exp(-0.3 * t), 1e-9);
  }
}

TEST(Cram, ZeroTimeIsIdentity) {
  auto A = linearChain({1.0, 0.5, 0.0});
  Eigen::VectorXd n0(3);
  n0 << 0.3, 0.6, 0.1;
  Eigen::VectorXd n = cramSolve(A, n0, 0.0, CramOrder::CRAM48);
  EXPECT_NEAR(n(0), 0.3, 1e-12);
  EXPECT_NEAR(n(1), 0.6, 1e-12);
  EXPECT_NEAR(n(2), 0.1, 1e-12);
}

TEST(Cram, StableSystemDoesNotChange) {
  Eigen::SparseMatrix<double> A(2, 2);  // all-zero matrix: nothing decays
  Eigen::VectorXd n0(2);
  n0 << 1.0, 2.0;
  Eigen::VectorXd n = cramSolve(A, n0, 1e6, CramOrder::CRAM48);
  EXPECT_NEAR(n(0), 1.0, 1e-12);
  EXPECT_NEAR(n(1), 2.0, 1e-12);
}

TEST(Cram, TwoStepChainMatchesBateman) {
  std::vector<double> lam{1.0e-2, 5.0e-3, 0.0};  // A -> B -> C(stable)
  auto A = linearChain(lam);
  double t = 60.0;
  auto ref = batemanLinearChain(lam, t);
  expectClose(cramSolve(A, unitFirst(3), t, CramOrder::CRAM48), ref, kRel48, kAbs);
  expectClose(cramSolve(A, unitFirst(3), t, CramOrder::CRAM16), ref, kRel16, kAbs);
}

TEST(Cram, FourStepChainMatchesBateman) {
  std::vector<double> lam{0.07, 0.013, 0.004, 0.0};
  auto A = linearChain(lam);
  double t = 50.0;
  auto ref = batemanLinearChain(lam, t);
  expectClose(cramSolve(A, unitFirst(4), t, CramOrder::CRAM48), ref, kRel48, kAbs);
  expectClose(cramSolve(A, unitFirst(4), t, CramOrder::CRAM16), ref, kRel16, kAbs);
}

TEST(Cram, BranchingMatchesAnalytic) {
  // A -> B (BR 0.7), A -> C (BR 0.3); B and C stable.
  const double la = 0.05, bAB = 0.7, bAC = 0.3, t = 30.0;
  std::vector<Eigen::Triplet<double>> tr{{0, 0, -la}, {1, 0, bAB * la}, {2, 0, bAC * la}};
  Eigen::SparseMatrix<double> A(3, 3);
  A.setFromTriplets(tr.begin(), tr.end());

  auto r = batemanTwoBranch(la, bAB, bAC, t);
  Eigen::VectorXd n = cramSolve(A, unitFirst(3), t, CramOrder::CRAM48);
  EXPECT_NEAR(n(0), r.a, kRel48 * r.a);
  EXPECT_NEAR(n(1), r.b, kRel48 * r.b);
  EXPECT_NEAR(n(2), r.c, kRel48 * r.c);
}

TEST(Cram, MassConservedForNonFissionChain) {
  auto A = linearChain({0.02, 0.011, 0.006, 0.0});
  Eigen::VectorXd n0 = unitFirst(4);
  for (double t : {1.0, 100.0, 1.0e4}) {
    Eigen::VectorXd n = cramSolve(A, n0, t, CramOrder::CRAM48);
    EXPECT_NEAR(n.sum(), 1.0, 1e-12) << "t=" << t;
  }
}

TEST(Cram, StiffSystemRemainsAccurate) {
  // Decay constants spanning ~8 orders of magnitude.
  std::vector<double> lam{1.0e-1, 1.0e-5, 1.0e-9, 0.0};
  auto A = linearChain(lam);
  double t = 1.0e5;
  auto ref = batemanLinearChain(lam, t);
  Eigen::VectorXd n = cramSolve(A, unitFirst(4), t, CramOrder::CRAM48);
  EXPECT_TRUE(n.allFinite());
  expectClose(n, ref, 1e-6, 1e-9);  // looser: stiff
  EXPECT_NEAR(n.sum(), 1.0, 1e-10);
}

TEST(Cram, InvalidShapeThrows) {
  Eigen::SparseMatrix<double> nonsquare(2, 3);
  Eigen::VectorXd n0(2);
  n0 << 1.0, 2.0;
  EXPECT_THROW(cramSolve(nonsquare, n0, 1.0, CramOrder::CRAM48), std::invalid_argument);

  Eigen::SparseMatrix<double> square(3, 3);
  Eigen::VectorXd mismatched(2);
  mismatched << 1.0, 2.0;
  EXPECT_THROW(cramSolve(square, mismatched, 1.0, CramOrder::CRAM48), std::invalid_argument);
}

TEST(Cram, Orders16And48Agree) {
  auto A = linearChain({0.03, 0.017, 0.009, 0.002, 0.0});
  double t = 80.0;
  Eigen::VectorXd n16 = cramSolve(A, unitFirst(5), t, CramOrder::CRAM16);
  Eigen::VectorXd n48 = cramSolve(A, unitFirst(5), t, CramOrder::CRAM48);
  for (int i = 0; i < n48.size(); ++i)
    if (std::abs(n48(i)) > 1e-8) {
      EXPECT_NEAR(n16(i), n48(i), 1e-4 * std::abs(n48(i))) << "index " << i;
    }
}
