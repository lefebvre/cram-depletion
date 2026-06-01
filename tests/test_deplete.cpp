#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "depletion_fixture.hpp"

using namespace cram;
using namespace cram_test;

TEST(Deplete, ReactionProductTopology) {
  const Zai u238{92, 238, 0};
  EXPECT_EQ(reactionProduct(u238, ReactionType::NGamma), (Zai{92, 239, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::N2n), (Zai{92, 237, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::N3n), (Zai{92, 236, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::N4n), (Zai{92, 235, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::NAlpha), (Zai{90, 235, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::NProton), (Zai{91, 238, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::Fission), (Zai{92, 238, 0}));  // unchanged
}

TEST(Deplete, SetReactionsReplacesAndZeroSigmaIsInert) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  sys.setConstantFlux(1.0e14);

  // First a real (n,gamma); then replace it with a zero-cross-section channel.
  sys.setReactions(kU235, {{ReactionType::NGamma, kU236, 9.0}});
  Eigen::SparseMatrix<double> withXs = sys.assemble(initialPinComposition(chain));
  const int i = chain.indexOf(kU235);
  EXPECT_LT(withXs.coeff(i, i), 0.0);  // U235 removed by capture

  sys.setReactions(kU235, {{ReactionType::NGamma, kU236, 0.0}});  // replaces prior set
  Eigen::SparseMatrix<double> noXs = sys.assemble(initialPinComposition(chain));
  EXPECT_EQ(noXs.coeff(i, i), 0.0);  // zero rate -> no removal term
}

TEST(Deplete, ConstantFluxMakesMatrixCompositionIndependent) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(2.5e14);

  Eigen::VectorXd n0 = initialPinComposition(chain);
  Eigen::VectorXd n1 = 0.5 * n0;

  Eigen::SparseMatrix<double> A0 = sys.assemble(n0);
  Eigen::SparseMatrix<double> A1 = sys.assemble(n1);
  EXPECT_EQ(sys.fluxFor(n0), sys.fluxFor(n1));
  EXPECT_TRUE(A0.isApprox(A1, 1e-15));
}

TEST(Deplete, ConstantPowerHoldsPowerAsCompositionEvolves) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);

  Eigen::VectorXd n0 = initialPinComposition(chain);
  sys.setConstantFlux(3.0e14);
  const double power = sys.powerFor(n0);
  sys.setConstantPower(power);

  // The set power is reproduced at BOL.
  EXPECT_NEAR(sys.powerFor(n0), power, power * 1e-12);
  EXPECT_NEAR(sys.fluxFor(n0), 3.0e14, 3.0e14 * 1e-9);

  auto integ = makeIntegrator(IntegratorKind::CECM, sys.matrixBuilder());
  auto res = deplete(*integ, n0, std::vector<double>(30, 10.0 * 86400.0));
  const Eigen::VectorXd& nEnd = res.n.back();

  EXPECT_LT(nEnd(chain.indexOf(kU235)), n0(chain.indexOf(kU235)));  // U235 burned
  EXPECT_GT(nEnd(chain.indexOf(kPu239)), 0.0);                      // Pu239 bred in
  // The flux floats (here it drops: bred Pu239 has a much larger fission cross
  // section than U235, so less flux is needed to hold the same power), but the
  // power must stay pinned at the target throughout.
  EXPECT_GT(std::abs(sys.fluxFor(nEnd) - sys.fluxFor(n0)) / sys.fluxFor(n0), 0.01);
  EXPECT_NEAR(sys.powerFor(nEnd), power, power * 1e-9);
}

TEST(Deplete, NoFissileGivesZeroFluxUnderConstantPower) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantPower(1.0);

  Eigen::VectorXd n(chain.size());
  n.setZero();
  n(chain.indexOf(kI135)) = 1.0;  // no fissile present
  EXPECT_EQ(sys.fluxFor(n), 0.0);

  // assemble then reduces to the decay matrix; I135 still decays.
  auto integ = makeIntegrator(IntegratorKind::Predictor, sys.matrixBuilder());
  auto res = deplete(*integ, n, std::vector<double>(5, 3600.0));
  EXPECT_LT(res.n.back()(chain.indexOf(kI135)), 1.0);
}

TEST(Deplete, TrajectoryShape) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(3.0e14);

  Eigen::VectorXd n0 = initialPinComposition(chain);
  auto integ = makeIntegrator(IntegratorKind::CELI, sys.matrixBuilder());
  std::vector<double> dts(12, 25.0 * 86400.0);
  auto res = deplete(*integ, n0, dts);

  ASSERT_EQ(res.n.size(), dts.size() + 1);
  ASSERT_EQ(res.time.size(), dts.size() + 1);
  EXPECT_EQ(res.time.front(), 0.0);
  EXPECT_NEAR(res.time.back(), 12 * 25.0 * 86400.0, 1e-3);
  // U235 monotonically decreasing under irradiation.
  const int iU235 = chain.indexOf(kU235);
  for (std::size_t k = 1; k < res.n.size(); ++k)
    EXPECT_LT(res.n[k](iU235), res.n[k - 1](iU235));
}
