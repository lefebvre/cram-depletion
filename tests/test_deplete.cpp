#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "depletion_pin_fixture.hpp"

using namespace cram;
using namespace cram_test;

TEST(Deplete, SetReactionsReplacesAndZeroSigmaIsInert) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  sys.setConstantFlux(1.0e14);

  // First a real (n,gamma); then replace it with a zero-cross-section channel.
  sys.setReactions(kU235, {capture(kU236, 9.0)});
  Eigen::SparseMatrix<double> withXs = sys.assemble(initialPinComposition(chain));
  const int i = chain.indexOf(kU235);
  EXPECT_LT(withXs.coeff(i, i), 0.0);  // U235 removed by capture

  sys.setReactions(kU235, {capture(kU236, 0.0)});  // replaces prior set
  Eigen::SparseMatrix<double> noXs = sys.assemble(initialPinComposition(chain));
  EXPECT_EQ(noXs.coeff(i, i), 0.0);  // zero rate -> no removal term
  ASSERT_EQ(sys.reactions().size(), 1u);
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
  EXPECT_EQ(sys.normalization(), DepletionSystem::Normalization::ConstantPower);

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

// A non-zero power target with no fissile material present is a contradiction
// the caller needs to hear about: returning zero flux would let a mis-keyed
// fission cross section make an entire irradiation silently vanish.
TEST(Deplete, NoFissileUnderNonZeroConstantPowerThrows) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantPower(1.0);

  Eigen::VectorXd n(chain.size());
  n.setZero();
  n(chain.indexOf(kI135)) = 1.0;  // no fissile present
  EXPECT_THROW(sys.fluxFor(n), std::domain_error);
  EXPECT_THROW(sys.assemble(n), std::domain_error);
}

// A zero power target is legitimately zero flux: assemble() reduces to the
// decay matrix and the march is a pure decay.
TEST(Deplete, ZeroConstantPowerIsPureDecay) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantPower(0.0);

  Eigen::VectorXd n(chain.size());
  n.setZero();
  n(chain.indexOf(kI135)) = 1.0;
  EXPECT_EQ(sys.fluxFor(n), 0.0);
  EXPECT_EQ(sys.powerFor(n), 0.0);

  auto integ = makeIntegrator(IntegratorKind::Predictor, sys.matrixBuilder());
  auto res = deplete(*integ, n, std::vector<double>(5, 3600.0));
  EXPECT_LT(res.n.back()(chain.indexOf(kI135)), 1.0);
  EXPECT_NEAR(res.n.back().sum(), 1.0, 1e-12);  // decay only: atoms conserved
}

// A fission channel without an energy release would contribute fissions but no
// power, making the constant-power normalization silently wrong. Rejected both
// when the power is set after the channel and when the channel is set after
// the power.
TEST(Deplete, FissionWithoutQIsRejectedUnderConstantPower) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  const ReactionXS noQ{.type = ReactionType::Fission, .target = std::nullopt, .sigma = 38.0};

  sys.setReactions(kU235, {noQ});  // fine under the default constant-flux mode
  EXPECT_THROW(sys.setConstantPower(1.0), std::invalid_argument);
  EXPECT_EQ(sys.normalization(), DepletionSystem::Normalization::ConstantFlux);

  sys.setReactions(kU235, {fission(38.0)});
  sys.setConstantPower(1.0);
  EXPECT_THROW(sys.setReactions(kPu239, {noQ}), std::invalid_argument);
  EXPECT_NO_THROW(sys.setConstantPower(0.0));  // zero power needs no Q
}

TEST(Deplete, NegativeCrossSectionIsRejected) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  EXPECT_THROW(sys.setReactions(kU235, {capture(kU236, -1.0)}), std::invalid_argument);
  EXPECT_TRUE(sys.reactions().empty());
}

// An untracked product (OpenMC's "Nothing" target) or a product the chain never
// registered still consumes the parent; only the production term is absent.
TEST(Deplete, UntrackedProductConsumesParentOnly) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  sys.setConstantFlux(1.0e14);
  const ReactionXS lost{.type = ReactionType::NAlpha, .target = std::nullopt, .sigma = 2.0};
  const ReactionXS unregistered{
      .type = ReactionType::N2n, .target = Zai{.z = 92, .a = 234}, .sigma = 3.0};
  sys.setReactions(kU235, {lost, unregistered});

  Eigen::SparseMatrix<double> A = sys.assemble(initialPinComposition(chain));
  const int i = chain.indexOf(kU235);
  const double expected = -(2.0 + 3.0) * 1e-24 * 1.0e14;
  EXPECT_NEAR(A.coeff(i, i), expected, 1e-12 * std::abs(expected));
  double colSum = 0.0;
  for (int r = 0; r < A.rows(); ++r)
    colSum += A.coeff(r, i);
  EXPECT_NEAR(colSum, A.coeff(i, i), 1e-12 * std::abs(expected)) << "no production anywhere";
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
