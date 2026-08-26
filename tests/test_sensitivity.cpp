#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cram/adjoint.hpp"
#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/deplete.hpp"
#include "depletion_pin_fixture.hpp"

using namespace cram;
using namespace cram_test;

namespace {

constexpr double kDay = 86400.0;

// A varying flux schedule with a shutdown, so the intervals genuinely differ.
const std::vector<double> kDts = {5.0 * kDay, 20.0 * kDay, 3.0 * kDay, 30.0 * kDay};
const std::vector<double> kFlux = {3.0e14, 1.0e14, 0.0, 2.0e14};

// Options fine enough for the finite-difference checks: under flux Xe-135 is
// removed in ~20 minutes, a tiny fraction of a 30-day interval, so the end
// pieces must be graded down to that scale.
const SensitivityOptions kFine{.gaussPoints = 5, .subIntervals = 16, .endRefinements = 12};

// The pin problem with one decay constant and one capture cross section
// exposed as parameters, so a perturbed copy can be built for finite
// differences.
struct Problem {
  double xeHalfLife = 32904.0;
  double u238Capture = 0.9;

  DepletionChain chain;
  std::vector<Eigen::SparseMatrix<double>> A;

  void build() {
    chain = buildPinChain();
    chain.setDecay(
        kXe135,
        DecayData{.halfLife = xeHalfLife,
                  .modes = {DecayMode{
                      .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
    DepletionSystem sys(chain);
    configurePinReactions(sys);
    sys.setReactions(kU238, {fission(0.05), capture(kU239, u238Capture)});
    sys.setConstantFlux(0.0);
    A = intervalMatrices(sys, kFlux);
  }

  // The system rebuilt on the current chain, for the contraction.
  DepletionSystem system() const {
    DepletionSystem sys(chain);
    configurePinReactions(sys);
    sys.setReactions(kU238, {fission(0.05), capture(kU239, u238Capture)});
    sys.setConstantFlux(0.0);
    return sys;
  }

  double finalResponse(const Eigen::VectorXd& w) const {
    return w.dot(depleteLinear(A, kDts, initialPinComposition(chain)).n.back());
  }
  double integratedResponse(const Eigen::VectorXd& w) const {
    return w.dot(integratedInventory(A, kDts, initialPinComposition(chain)));
  }
};

Eigen::VectorXd unit(const DepletionChain& chain, const Zai& z) {
  Eigen::VectorXd w = Eigen::VectorXd::Zero(chain.size());
  w(chain.indexOf(z)) = 1.0;
  return w;
}

}  // namespace

// dn/dt = A n + s for one decaying nuclide with a constant feed has the closed
// form n(t) = n0 e^{-lambda t} + s/lambda (1 - e^{-lambda t}).
TEST(Sensitivity, SourceGeneratorSolvesTheInhomogeneousEquation) {
  const double lambda = 2.0e-5, s0 = 3.0, n0 = 1.0, t = 5.0e4;
  Eigen::SparseMatrix<double> A(1, 1);
  A.insert(0, 0) = -lambda;
  A.makeCompressed();
  Eigen::VectorXd s(1);
  s << s0;
  Eigen::VectorXd x(2);
  x << n0, 1.0;
  const Eigen::VectorXd y = cramSolve(sourceGenerator(A, s), x, t);
  const double e = std::exp(-lambda * t);
  EXPECT_NEAR(y(0), n0 * e + s0 / lambda * (1.0 - e), 1e-12 * (s0 / lambda));
  EXPECT_NEAR(y(1), 1.0, 1e-13);
}

// The augmented generator's integral against the closed-form Bateman integrals
// of I-135 -> Xe-135 -> Cs-135, accumulated across two intervals:
//   int N_I  = (1 - e1) / l1
//   int N_Xe = l1/(l2 - l1) * ((1 - e1)/l1 - (1 - e2)/l2)
//   int N_Cs = T - int N_I - int N_Xe
TEST(Sensitivity, IntegratedInventoryMatchesAnalyticBateman) {
  DepletionChain chain;
  const Zai I{53, 135, 0}, Xe{54, 135, 0}, Cs{55, 135, 0};
  const auto beta = [](double halfLife) {
    return DecayData{
        .halfLife = halfLife,
        .modes = {DecayMode{.rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}};
  };
  chain.setDecay(I, beta(23652.0));
  chain.setDecay(Xe, beta(32904.0));
  chain.close();
  const double l1 = chain.decay(I)->decayConstant, l2 = chain.decay(Xe)->decayConstant;
  const Eigen::SparseMatrix<double> A = chain.decayMatrix();

  const double T = 40.0 * 3600.0;
  const std::vector<double> dts = {0.3 * T, 0.7 * T};
  const std::vector<Eigen::SparseMatrix<double>> As = {A, A};
  Eigen::VectorXd n0 = Eigen::VectorXd::Zero(chain.size());
  n0(chain.indexOf(I)) = 1.0;
  const Eigen::VectorXd integral = integratedInventory(As, dts, n0);

  const double e1 = std::exp(-l1 * T), e2 = std::exp(-l2 * T);
  const double intI = (1.0 - e1) / l1;
  const double intXe = l1 / (l2 - l1) * ((1.0 - e1) / l1 - (1.0 - e2) / l2);
  EXPECT_NEAR(integral(chain.indexOf(I)), intI, 1e-12 * intI);
  EXPECT_NEAR(integral(chain.indexOf(Xe)), intXe, 1e-12 * intXe);
  EXPECT_NEAR(integral(chain.indexOf(Cs)), T - intI - intXe, 1e-12 * T);
}

// The integrated adjoint's initial value is the exact gradient of the
// time-integrated response with respect to n0.
TEST(Sensitivity, IntegratedAdjointIsTheGradientOfTheIntegratedResponse) {
  Problem p;
  p.build();
  const Eigen::VectorXd w = unit(p.chain, kXe135) + 0.5 * unit(p.chain, kCs135);
  const AdjointResult adj = adjointDepleteIntegrated(p.A, kDts, w);
  ASSERT_TRUE(adj.integratedWeight.has_value());
  EXPECT_EQ(adj.nStar.back().norm(), 0.0);

  const Eigen::VectorXd n0 = initialPinComposition(p.chain);
  const double R = w.dot(integratedInventory(p.A, kDts, n0));
  EXPECT_NEAR(adj.nStar.front().dot(n0), R, 1e-9 * R) << "duality";

  Eigen::VectorXd delta = Eigen::VectorXd::Zero(n0.size());
  delta(p.chain.indexOf(kPu239)) = 1.0e-5;
  const double Rp = w.dot(integratedInventory(p.A, kDts, n0 + delta));
  EXPECT_NEAR(Rp - R, adj.nStar.front().dot(delta), 1e-9 * std::abs(Rp - R));
}

// dR/d(lambda_Xe) from the adjoint sensitivities against a central finite
// difference of the forward solution.
TEST(Sensitivity, DecayConstantSensitivityMatchesFiniteDifference) {
  Problem p;
  p.build();
  const Eigen::VectorXd w = unit(p.chain, kCs135);  // Cs-135 at end of schedule
  const int iXe = p.chain.indexOf(kXe135);

  const DepletionResult fwd = depleteLinear(p.A, kDts, initialPinComposition(p.chain));
  const AdjointResult adj = adjointDeplete(p.A, kDts, w);
  const auto S = rateSensitivities(fwd, adj, p.A, kDts, kFine);
  const Eigen::VectorXd dRdLambda = decayConstantSensitivities(p.chain, S);

  const double lambda = kLn2 / p.xeHalfLife;
  const double eps = 1.0e-5;
  Problem plus = p, minus = p;
  plus.xeHalfLife = kLn2 / (lambda * (1.0 + eps));
  minus.xeHalfLife = kLn2 / (lambda * (1.0 - eps));
  plus.build();
  minus.build();
  const double fd = (plus.finalResponse(w) - minus.finalResponse(w)) / (2.0 * eps * lambda);

  EXPECT_NE(fd, 0.0);
  EXPECT_NEAR(dRdLambda(iXe), fd, 1e-6 * std::abs(fd));
  // Nuclides without decay data carry no sensitivity.
  EXPECT_EQ(dRdLambda(p.chain.indexOf(kU235)), 0.0);
}

// dR/d(sigma) for U-238 capture against a central finite difference, for a
// final-time response and for a time-integrated one.
TEST(Sensitivity, CaptureCrossSectionSensitivityMatchesFiniteDifference) {
  Problem p;
  p.build();
  const Eigen::VectorXd w = unit(p.chain, kPu239);
  const DepletionResult fwd = depleteLinear(p.A, kDts, initialPinComposition(p.chain));

  const double eps = 1.0e-5;
  Problem plus = p, minus = p;
  plus.u238Capture = p.u238Capture * (1.0 + eps);
  minus.u238Capture = p.u238Capture * (1.0 - eps);
  plus.build();
  minus.build();
  const double dSigma = 2.0 * eps * p.u238Capture;

  auto captureEntry = [&](const std::vector<ReactionSensitivity>& rs) {
    for (const auto& r : rs)
      if (r.parent == kU238 && r.type == ReactionType::NGamma)
        return r.dRdSigma;
    ADD_FAILURE() << "no U-238 capture entry";
    return 0.0;
  };

  {
    const AdjointResult adj = adjointDeplete(p.A, kDts, w);
    const auto S = rateSensitivities(fwd, adj, p.A, kDts, kFine);
    const double got = captureEntry(reactionSensitivities(p.system(), S, kFlux));
    const double fd = (plus.finalResponse(w) - minus.finalResponse(w)) / dSigma;
    EXPECT_NE(fd, 0.0);
    EXPECT_NEAR(got, fd, 1e-6 * std::abs(fd)) << "final-time response";
  }
  {
    const AdjointResult adj = adjointDepleteIntegrated(p.A, kDts, w);
    const auto S = rateSensitivities(fwd, adj, p.A, kDts, kFine);
    const double got = captureEntry(reactionSensitivities(p.system(), S, kFlux));
    const double fd = (plus.integratedResponse(w) - minus.integratedResponse(w)) / dSigma;
    EXPECT_NE(fd, 0.0);
    EXPECT_NEAR(got, fd, 1e-6 * std::abs(fd)) << "time-integrated response";
  }
}

// Refining the quadrature converges: the fine rule and a finer one agree far
// more closely than the coarse rule agrees with either.
TEST(Sensitivity, QuadratureConverges) {
  Problem p;
  p.build();
  const Eigen::VectorXd w = unit(p.chain, kCs135);
  const DepletionResult fwd = depleteLinear(p.A, kDts, initialPinComposition(p.chain));
  const AdjointResult adj = adjointDeplete(p.A, kDts, w);
  const int iXe = p.chain.indexOf(kXe135);

  auto xeSensitivity = [&](SensitivityOptions o) {
    return decayConstantSensitivities(p.chain, rateSensitivities(fwd, adj, p.A, kDts, o))(iXe);
  };
  const double coarse = xeSensitivity({.gaussPoints = 3, .subIntervals = 1, .endRefinements = 0});
  const double defaults = xeSensitivity({});
  const double fine = xeSensitivity(kFine);
  const double finer = xeSensitivity({.gaussPoints = 5, .subIntervals = 32, .endRefinements = 14});
  EXPECT_NEAR(fine, finer, 1e-8 * std::abs(finer));
  EXPECT_GT(std::abs(coarse - finer), 100.0 * std::abs(fine - finer));
  // The defaults sit between: usable, and visibly improvable.
  EXPECT_LT(std::abs(defaults - finer), std::abs(coarse - finer));
}

TEST(Sensitivity, RejectsBadOptionsAndMismatchedInputs) {
  Problem p;
  p.build();
  const Eigen::VectorXd w = unit(p.chain, kCs135);
  const DepletionResult fwd = depleteLinear(p.A, kDts, initialPinComposition(p.chain));
  const AdjointResult adj = adjointDeplete(p.A, kDts, w);
  EXPECT_THROW(rateSensitivities(fwd, adj, p.A, kDts, {.gaussPoints = 6}), std::invalid_argument);
  EXPECT_THROW(rateSensitivities(fwd, adj, p.A, kDts, {.subIntervals = 0}), std::invalid_argument);
  EXPECT_THROW(rateSensitivities(fwd, adj, p.A, kDts, {.endRefinements = -1}),
               std::invalid_argument);
  const std::vector<double> shortDts(kDts.begin(), kDts.end() - 1);
  const std::vector<Eigen::SparseMatrix<double>> shortA(p.A.begin(), p.A.end() - 1);
  EXPECT_THROW(rateSensitivities(fwd, adj, shortA, shortDts), std::invalid_argument);
  const auto S = rateSensitivities(fwd, adj, p.A, kDts);
  EXPECT_THROW(reactionSensitivities(p.system(), S, {1.0}), std::invalid_argument);
  EXPECT_THROW(sourceGenerator(p.A[0], Eigen::VectorXd::Zero(2)), std::invalid_argument);
}
