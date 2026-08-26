#include <gtest/gtest.h>

#include <cmath>
#include <limits>
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

// A fissionable nuclide the chain carries no yield table for (a chain trimmed
// of its yields, or one whose delegated yields never resolved) must still be
// consumed by its fission channel. Dropping the channel whole leaves the
// parent untouched while fissionPowerWeight() keeps crediting its energy
// release, so a constant-power calculation holds a power the matrix never
// delivers and the material under-burns.
TEST(Deplete, FissionWithoutYieldsStillConsumesTheParent) {
  DepletionChain chain = buildPinChain();  // U238 gets a fission xs, but no yields
  ASSERT_EQ(chain.nearestYields(kU238, 0.0253), nullptr);
  DepletionSystem sys(chain);
  sys.setConstantFlux(1.0e14);
  sys.setReactions(kU238, {fission(0.05)});

  const Eigen::SparseMatrix<double> A = sys.assemble(initialPinComposition(chain));
  const int i = chain.indexOf(kU238);
  const double expected = -0.05 * 1e-24 * 1.0e14;
  EXPECT_NEAR(A.coeff(i, i), expected, 1e-12 * std::abs(expected));
  double colSum = 0.0;
  for (int r = 0; r < A.rows(); ++r)
    colSum += A.coeff(r, i);
  EXPECT_NEAR(colSum, A.coeff(i, i), 1e-12 * std::abs(expected)) << "no products tracked";
}

// Every fission the power normalization pays for is a fission the matrix
// performs: with only fission channels active, each parent's diagonal is its
// fission rate, and those rates must reproduce the requested power exactly.
// U238 here has no yield table, which is the case that used to go missing.
TEST(Deplete, ConstantPowerMatchesTheFissionsTheMatrixPerforms) {
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  sys.setReactions(kU235, {fission(38.0)});
  sys.setReactions(kU238, {fission(0.05)});
  const Eigen::VectorXd n = initialPinComposition(chain);
  sys.setConstantFlux(1.0e14);
  const double target = sys.powerFor(n);
  sys.setConstantPower(target);

  const Eigen::SparseMatrix<double> A = sys.assemble(n);
  constexpr double kEvToJoule = 1.602176634e-19;
  double delivered = 0.0;
  for (const Zai& z : {kU235, kU238}) {
    const int i = chain.indexOf(z);
    delivered += n(i) * -A.coeff(i, i) * kFissionQ * kEvToJoule;
  }
  EXPECT_NEAR(delivered, target, 1e-12 * target);
}

// A negative flux or power runs every reaction backwards; a non-finite one
// poisons the factorization of every matrix assembled from it. Neither has a
// reading, and both are rejected without disturbing the current mode.
TEST(Deplete, NonFiniteOrNegativeNormalizationIsRejected) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(2.0e14);

  for (double bad : {-1.0, nan, inf}) {
    EXPECT_THROW(sys.setConstantFlux(bad), std::invalid_argument);
    EXPECT_THROW(sys.setConstantPower(bad), std::invalid_argument);
  }
  EXPECT_EQ(sys.normalization(), DepletionSystem::Normalization::ConstantFlux);
  EXPECT_DOUBLE_EQ(sys.fluxFor(initialPinComposition(chain)), 2.0e14);
}

// The decay half of the matrix is cached at construction, the reaction half is
// sized from the chain as it stands at assemble(). Growing the chain in
// between -- which is exactly what the constructor's own "call close() before
// building the matrix" warning invites -- must not mix the two.
TEST(Deplete, ChainModifiedAfterConstructionIsRefused) {
  const Zai kCm242{.z = 96, .a = 242, .i = 0};  // alpha -> Pu238, absent here
  DepletionChain chain = buildPinChain();
  DepletionSystem sys(chain);
  configurePinReactions(sys);
  sys.setConstantFlux(1.0e14);
  const Eigen::VectorXd n = initialPinComposition(chain);
  EXPECT_NO_THROW(sys.assemble(n));

  chain.setDecay(
      kCm242, DecayData{.halfLife = 1.4e7,
                        .modes = {DecayMode{
                            .rtyp = 4.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  EXPECT_THROW(sys.assemble(n), std::logic_error);

  EXPECT_EQ(chain.close(), 1);  // registers Pu238
  EXPECT_THROW(sys.assemble(n), std::logic_error);

  sys.refreshChain();
  Eigen::VectorXd grown = Eigen::VectorXd::Zero(chain.size());
  grown.head(n.size()) = n;
  const Eigen::SparseMatrix<double> A = sys.assemble(grown);
  EXPECT_EQ(A.rows(), chain.size());
  const int cm = chain.indexOf(kCm242);
  EXPECT_LT(A.coeff(cm, cm), 0.0);
  EXPECT_GT(A.coeff(chain.indexOf({94, 238, 0}), cm), 0.0);
}
