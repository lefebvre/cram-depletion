#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cram/cram.hpp"
#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "depletion_fixture.hpp"

using namespace cram;
using namespace cram_test;

namespace {

std::vector<double> uniformSchedule(double total, int steps) {
  return std::vector<double>(static_cast<std::size_t>(steps), total / steps);
}

constexpr double kDay = 86400.0;

const std::vector<IntegratorKind> kAllKinds = {IntegratorKind::Predictor, IntegratorKind::CECM,
                                               IntegratorKind::CELI, IntegratorKind::LEQI,
                                               IntegratorKind::CF4};

}  // namespace

// When the flux (and hence A) is composition-independent, every scheme must
// collapse to the exact exp(A*T) and therefore agree with a single CRAM solve.
TEST(Integrator, ConstantFluxAllSchemesAreExact) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(3.0e14);

  Eigen::VectorXd n0 = initialPinComposition(chain);
  const double T = 120.0 * kDay;

  // A is constant, so exp(A*T) n0 is the exact answer.
  Eigen::SparseMatrix<double> A = sys.assemble(n0);
  Eigen::VectorXd exact = cramSolve(A, n0, T, CramOrder::CRAM48);

  for (IntegratorKind kind : kAllKinds) {
    auto integ = makeIntegrator(kind, sys.matrixBuilder());
    auto res = deplete(*integ, n0, uniformSchedule(T, 8));
    const Eigen::VectorXd& got = res.n.back();
    for (int i = 0; i < got.size(); ++i) {
      const double tol = 1e-9 * std::max(1.0, std::abs(exact(i)));
      EXPECT_NEAR(got(i), exact(i), tol) << integratorName(kind) << " index " << i;
    }
  }
}

// Pure decay (zero flux) of a beta-chain conserves the total atom count, and
// the marched solution must match a single exp(A*T).
TEST(Integrator, DecayOnlyConservesAtomsAndMatchesDirect) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(0.0);  // no reactions -> decay only

  Eigen::VectorXd n0(chain.size());
  n0.setZero();
  n0(chain.indexOf(kI135)) = 1.0;  // I135 -> Xe135 -> Cs135 (stable)
  const double T = 20.0 * 3600.0;

  Eigen::SparseMatrix<double> A = sys.assemble(n0);
  Eigen::VectorXd direct = cramSolve(A, n0, T, CramOrder::CRAM48);

  auto integ = makeIntegrator(IntegratorKind::CECM, sys.matrixBuilder());
  auto res = deplete(*integ, n0, uniformSchedule(T, 10));

  EXPECT_NEAR(res.n.back().sum(), 1.0, 1e-9);  // mass conservation
  for (int i = 0; i < direct.size(); ++i)
    EXPECT_NEAR(res.n.back()(i), direct(i), 1e-9) << "index " << i;
}

// Observed time-step convergence order, mirroring the paper's section 4.5 /
// Figs 15-18: predictor is first order, CE/CM, CE/LI, LE/QI are second order,
// CF4 is high order. Error is measured against a finely-stepped CF4 reference.
TEST(Integrator, ConvergenceOrderMatchesScheme) {
  DepletionChain chain = buildPinChain();
  const Eigen::VectorXd n0 = initialPinComposition(chain);
  const double T = 150.0 * kDay;
  const int iU235 = chain.indexOf(kU235);

  // Constant-power normalization: hold the power that gives ~1e14 n/cm^2/s at
  // BOL. The flux then floats as the composition evolves, making A genuinely
  // composition-dependent so the integrator order is exercised. The flux is
  // kept modest so the schemes reach their asymptotic regime at these step
  // counts (a much harder burn shows non-monotonic pre-asymptotic error).
  DepletionSystem psys(chain);
  configurePinReactions(psys);
  psys.setConstantFlux(1.0e14);
  psys.setConstantPower(psys.powerFor(n0));

  // High-accuracy CF4 reference (4th order, well below the errors we measure).
  auto refInteg = makeIntegrator(IntegratorKind::CF4, psys.matrixBuilder());
  const double truth = deplete(*refInteg, n0, uniformSchedule(T, 320)).n.back()(iU235);

  const std::vector<int> stepCounts = {16, 32, 64};  // s0, s1, s2

  auto errorsFor = [&](IntegratorKind kind) {
    std::vector<double> errs;
    for (int s : stepCounts) {
      auto integ = makeIntegrator(kind, psys.matrixBuilder());
      double x = deplete(*integ, n0, uniformSchedule(T, s)).n.back()(iU235);
      errs.push_back(std::abs(x - truth) / std::abs(truth));
    }
    return errs;
  };
  auto order = [](double eCoarse, double eFine) { return std::log2(eCoarse / eFine); };

  const auto ePred = errorsFor(IntegratorKind::Predictor);
  const auto eCecm = errorsFor(IntegratorKind::CECM);
  const auto eCeli = errorsFor(IntegratorKind::CELI);
  const auto eLeqi = errorsFor(IntegratorKind::LEQI);
  const auto eCf4 = errorsFor(IntegratorKind::CF4);

  // Predictor is first order; the finest pair sits in the asymptotic regime.
  EXPECT_NEAR(order(ePred[1], ePred[2]), 1.0, 0.2);

  // CE/CM, CE/LI, LE/QI are (at least) second order.
  EXPECT_GT(order(eCecm[1], eCecm[2]), 1.7) << "cecm";
  EXPECT_GT(order(eCeli[1], eCeli[2]), 1.7) << "celi";
  EXPECT_GT(order(eLeqi[1], eLeqi[2]), 1.7) << "leqi";

  // CF4 is high order (measured on the coarser pair, above the roundoff floor).
  EXPECT_GT(order(eCf4[0], eCf4[1]), 3.0) << "cf4";

  // Accuracy hierarchy at the finest step, mirroring the paper's finding that
  // the predictor is far less accurate than the higher-order schemes.
  EXPECT_LT(eCecm[2] * 50.0, ePred[2]) << "cecm vs predictor";
  EXPECT_LT(eCeli[2] * 50.0, ePred[2]) << "celi vs predictor";
  EXPECT_LT(eLeqi[2] * 50.0, ePred[2]) << "leqi vs predictor";
  EXPECT_LT(eCf4[2], eCecm[2]) << "cf4 most accurate";
}

TEST(Integrator, NameRoundTrip) {
  for (IntegratorKind kind : kAllKinds)
    EXPECT_EQ(integratorKindFromName(integratorName(kind)), kind);
  EXPECT_EQ(integratorKindFromName("CE/CM"), IntegratorKind::CECM);
  EXPECT_EQ(integratorKindFromName("LE-QI"), IntegratorKind::LEQI);
  EXPECT_THROW(integratorKindFromName("bogus"), std::invalid_argument);
}
