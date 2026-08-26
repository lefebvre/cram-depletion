#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cram/adjoint.hpp"
#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "cyclic_chain.hpp"
#include "depletion_pin_fixture.hpp"

using namespace cram;
using namespace cram_test;

namespace {

constexpr double kDay = 86400.0;

// I-135 -> Xe-135 -> Cs-135 (stable), with the fixture's half-lives.
struct IodineChain {
  DepletionChain chain;
  int iI, iXe, iCs;
  double lamI, lamXe;
  Eigen::SparseMatrix<double> A;

  IodineChain() {
    const Zai I{53, 135, 0}, Xe{54, 135, 0}, Cs{55, 135, 0};
    chain.setDecay(
        I, DecayData{.halfLife = 23652.0,
                     .modes = {DecayMode{
                         .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
    chain.setDecay(
        Xe, DecayData{.halfLife = 32904.0,
                      .modes = {DecayMode{
                          .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
    chain.close();
    iI = chain.indexOf(I);
    iXe = chain.indexOf(Xe);
    iCs = chain.indexOf(Cs);
    lamI = chain.decay(I)->decayConstant;
    lamXe = chain.decay(Xe)->decayConstant;
    A = chain.decayMatrix();
  }
};

}  // namespace

TEST(Adjoint, TransposedIsCompressedAndSwapsIndices) {
  const CyclicChain c = buildCyclicChain(40);
  const Eigen::SparseMatrix<double> At = transposed(c.A);
  EXPECT_TRUE(At.isCompressed());
  EXPECT_EQ(At.nonZeros(), c.A.nonZeros());
  for (int k = 0; k < c.A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(c.A, k); it; ++it)
      EXPECT_EQ(At.coeff(it.col(), it.row()), it.value());
}

// <w, exp(A dt) n0> == <exp(A^T dt) w, n0> on a cyclic burnup matrix: the
// adjoint solve is the transpose of the forward one.
TEST(Adjoint, DualityIdentityOnCyclicChain) {
  const CyclicChain c = buildCyclicChain(60);
  Eigen::VectorXd w(c.A.rows());
  for (int i = 0; i < w.size(); ++i)
    w(i) = 0.5 + 0.01 * i;  // an arbitrary, dense response weight
  const double dt = 3.0 * kDay;

  const Eigen::VectorXd n = cramSolve(c.A, c.n0, dt, CramOrder::CRAM48);
  const Eigen::VectorXd nStar = cramSolveAdjoint(c.A, w, dt, CramOrder::CRAM48);
  const double forward = w.dot(n);
  const double adjoint = nStar.dot(c.n0);
  EXPECT_NEAR(adjoint, forward, 1e-12 * std::abs(forward));
}

// Closed-form importance for the Cs-135 inventory at T: the probability that
// an atom present at t = 0 has become Cs-135 by T.
//   Xe-135: 1 - exp(-lamXe T)
//   I-135:  1 - exp(-lamI T) - lamI/(lamXe - lamI) (exp(-lamI T) - exp(-lamXe T))
//   Cs-135: 1
TEST(Adjoint, MatchesAnalyticBatemanImportance) {
  const IodineChain c;
  const double T = 30.0 * 3600.0;
  Eigen::VectorXd w = Eigen::VectorXd::Zero(c.chain.size());
  w(c.iCs) = 1.0;

  const Eigen::VectorXd nStar = cramSolveAdjoint(c.A, w, T);
  const double eI = std::exp(-c.lamI * T), eXe = std::exp(-c.lamXe * T);
  EXPECT_NEAR(nStar(c.iXe), 1.0 - eXe, 1e-13);
  EXPECT_NEAR(nStar(c.iI), 1.0 - eI - c.lamI / (c.lamXe - c.lamI) * (eI - eXe), 1e-13);
  EXPECT_NEAR(nStar(c.iCs), 1.0, 1e-13);
}

// dR/dn0 from the adjoint equals the change in R from perturbing n0, exactly,
// because the system is linear.
TEST(Adjoint, InitialImportanceIsTheGradientOfTheResponse) {
  const IodineChain c;
  const double T = 12.0 * 3600.0;
  Eigen::VectorXd w = Eigen::VectorXd::Zero(c.chain.size());
  w(c.iXe) = 1.0;  // response: Xe-135 inventory at T

  Eigen::VectorXd n0 = Eigen::VectorXd::Zero(c.chain.size());
  n0(c.iI) = 2.0;
  n0(c.iXe) = 0.5;
  const double R = w.dot(cramSolve(c.A, n0, T));

  const Eigen::VectorXd grad = cramSolveAdjoint(c.A, w, T);
  Eigen::VectorXd delta = Eigen::VectorXd::Zero(c.chain.size());
  delta(c.iI) = 0.3;
  const double Rp = w.dot(cramSolve(c.A, n0 + delta, T));
  EXPECT_NEAR(Rp - R, grad.dot(delta), 1e-12);
}

// A piecewise-constant schedule with equal flux on every interval is one
// interval in disguise: the march must agree with the single solve, forward
// and backward, and the duality identity must hold across the whole march.
TEST(Adjoint, PiecewiseMarchMatchesSingleIntervalUnderEqualFlux) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(2.0e14);

  const std::vector<double> dts(6, 20.0 * kDay);
  const std::vector<double> flux(6, 2.0e14);
  const auto A = intervalMatrices(sys, flux);
  ASSERT_EQ(A.size(), 6u);
  const double T = 6 * 20.0 * kDay;

  const Eigen::VectorXd n0 = initialPinComposition(chain);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(chain.size());
  w(chain.indexOf(kPu239)) = 1.0;

  const DepletionResult fwd = depleteLinear(A, dts, n0);
  ASSERT_EQ(fwd.n.size(), 7u);
  EXPECT_NEAR(fwd.time.back(), T, 1e-6);
  const Eigen::VectorXd direct = cramSolve(A[0], n0, T);
  for (int i = 0; i < direct.size(); ++i)
    EXPECT_NEAR(fwd.n.back()(i), direct(i), 1e-9 * std::max(1.0, std::abs(direct(i)))) << i;

  const AdjointResult adj = adjointDeplete(A, dts, w);
  ASSERT_EQ(adj.nStar.size(), 7u);
  EXPECT_EQ(adj.time, fwd.time);
  EXPECT_TRUE(adj.nStar.back().isApprox(w));
  const Eigen::VectorXd directStar = cramSolveAdjoint(A[0], w, T);
  for (int i = 0; i < directStar.size(); ++i)
    EXPECT_NEAR(adj.nStar.front()(i), directStar(i), 1e-9 * std::max(1.0, std::abs(directStar(i))))
        << i;

  // <n*(t_k), n(t_k)> is invariant along the march.
  const double R = w.dot(fwd.n.back());
  for (std::size_t k = 0; k < fwd.n.size(); ++k)
    EXPECT_NEAR(adj.nStar[k].dot(fwd.n[k]), R, 1e-9 * R) << "point " << k;
}

// With a genuinely varying flux the march is not a single exponential, but
// duality still holds and the linear march agrees with the predictor
// integrator driven by the same matrices.
TEST(Adjoint, VaryingFluxKeepsDualityAndMatchesPredictor) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(0.0);

  const std::vector<double> dts = {5.0 * kDay, 30.0 * kDay, 2.0 * kDay, 60.0 * kDay};
  const std::vector<double> flux = {3.0e14, 1.0e14, 0.0, 2.5e14};  // includes a shutdown
  const auto A = intervalMatrices(sys, flux);

  const Eigen::VectorXd n0 = initialPinComposition(chain);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(chain.size());
  w(chain.indexOf(kXe135)) = 1.0;
  w(chain.indexOf(kCs135)) = 0.25;

  const DepletionResult fwd = depleteLinear(A, dts, n0);
  const AdjointResult adj = adjointDeplete(A, dts, w);
  const double R = w.dot(fwd.n.back());
  EXPECT_NEAR(adj.nStar.front().dot(n0), R, 1e-9 * R);

  // The same schedule through the generic integrator, one interval at a time.
  std::size_t k = 0;
  auto integ =
      makeIntegrator(IntegratorKind::Predictor, [&](const Eigen::VectorXd&) { return A[k]; });
  Eigen::VectorXd n = n0;
  for (k = 0; k < dts.size(); ++k)
    n = integ->step(n, dts[k]);
  for (int i = 0; i < n.size(); ++i)
    EXPECT_NEAR(fwd.n.back()(i), n(i), 1e-9 * std::max(1.0, std::abs(n(i)))) << i;
}

TEST(Adjoint, IntervalMatricesRefuseConstantPower) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantPower(1.0);
  EXPECT_THROW(intervalMatrices(sys, {1.0e14}), std::invalid_argument);
}

TEST(Adjoint, MismatchedScheduleIsRejected) {
  const CyclicChain c = buildCyclicChain(10);
  const std::vector<Eigen::SparseMatrix<double>> A = {c.A, c.A};
  EXPECT_THROW(depleteLinear(A, {1.0}, c.n0), std::invalid_argument);
  EXPECT_THROW(adjointDeplete(A, {1.0, 2.0, 3.0}, c.n0), std::invalid_argument);
}
