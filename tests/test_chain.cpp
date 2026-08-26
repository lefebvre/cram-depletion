#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "cram/chain.hpp"
#include "cram/endf_reader.hpp"

using namespace cram;

// ---------------------------------------------------------------------------
// Registration / lookup
// ---------------------------------------------------------------------------
TEST(Chain, AddIsIdempotent) {
  DepletionChain c;
  int i = c.add({54, 135, 0});
  int j = c.add({54, 135, 0});
  EXPECT_EQ(i, j);
  EXPECT_EQ(c.size(), 1);
  c.add({55, 135, 0});
  EXPECT_EQ(c.size(), 2);
}

TEST(Chain, IndexOfAbsent) {
  DepletionChain c;
  c.add({1, 1, 0});
  EXPECT_EQ(c.indexOf({1, 1, 0}), 0);
  EXPECT_EQ(c.indexOf({2, 4, 0}), -1);
}

// ---------------------------------------------------------------------------
// Deterministic ordering
// ---------------------------------------------------------------------------
namespace {

// Five beta- emitters whose daughters are deliberately left unregistered, so
// close() is what introduces them.
const std::vector<Zai> kParents = {
    {40, 97, 0}, {52, 133, 0}, {36, 89, 0}, {58, 143, 0}, {44, 105, 0},
};

DecayData betaMinus(double halfLife) {
  return DecayData{
      .halfLife = halfLife,
      .modes = {DecayMode{.rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}};
}

// The nuclides close() appended, in the order it appended them.
std::vector<std::int64_t> daughtersAddedByClose(const std::vector<Zai>& insertionOrder) {
  DepletionChain c;
  for (const Zai& z : insertionOrder)
    c.setDecay(z, betaMinus(1000.0));
  const int before = c.size();
  c.close();

  std::vector<std::int64_t> added;
  for (int i = before; i < c.size(); ++i)
    added.push_back(c.nuclides()[static_cast<std::size_t>(i)].key());
  return added;
}

}  // namespace

// close() walks the decaying nuclides to find daughters to register, and add()
// hands out matrix indices in call order. Walking an unordered_map would let the
// hash function decide that order -- and therefore decide the matrix layout, the
// LU pivoting, and the last bits of every result -- differently on different
// standard libraries. Feeding the same nuclides in a different sequence must not
// change the order close() introduces their daughters in.
TEST(Chain, CloseAddsDaughtersInAnInsertionOrderIndependentOrder) {
  std::vector<Zai> shuffled = kParents;
  std::rotate(shuffled.begin(), shuffled.begin() + 2, shuffled.end());
  std::vector<Zai> reversed(kParents.rbegin(), kParents.rend());

  const std::vector<std::int64_t> a = daughtersAddedByClose(kParents);
  const std::vector<std::int64_t> b = daughtersAddedByClose(shuffled);
  const std::vector<std::int64_t> c = daughtersAddedByClose(reversed);

  ASSERT_EQ(a.size(), kParents.size()) << "each parent should contribute one daughter";
  EXPECT_EQ(a, b);
  EXPECT_EQ(a, c);
  // Ascending ZAI, which is what sorting the parent keys before the walk gives.
  EXPECT_TRUE(std::is_sorted(a.begin(), a.end()));
}

// The triplet stream reaches setFromTriplets(), which sums duplicate entries in
// the order it receives them. Hash-ordered iteration would make that summation
// order, and so the floating-point result, implementation-defined.
TEST(Chain, DecayMatrixIsIdenticalRegardlessOfInsertionOrder) {
  auto build = [](const std::vector<Zai>& order) {
    DepletionChain c;
    for (const Zai& z : order)
      c.setDecay(z, betaMinus(1000.0));
    c.close();
    return c.decayMatrix();
  };

  std::vector<Zai> reversed(kParents.rbegin(), kParents.rend());
  const Eigen::SparseMatrix<double> forward = build(kParents);
  const Eigen::SparseMatrix<double> backward = build(reversed);

  // Parents keep their own insertion indices, so the matrices are not expected
  // to be identical -- but they must have the same shape and nonzero count, and
  // the trace is order-independent, so it should agree to the last bit.
  ASSERT_EQ(forward.rows(), backward.rows());
  EXPECT_EQ(forward.nonZeros(), backward.nonZeros());
  double traceF = 0.0;
  double traceB = 0.0;
  for (int i = 0; i < forward.rows(); ++i) {
    traceF += forward.coeff(i, i);
    traceB += backward.coeff(i, i);
  }
  EXPECT_DOUBLE_EQ(traceF, traceB);
}

// ---------------------------------------------------------------------------
// Decay data
// ---------------------------------------------------------------------------
TEST(Chain, DecayConstantFromHalfLife) {
  DepletionChain c;
  c.setDecay({54, 135, 0}, DecayData{.halfLife = 100.0});
  const DecayData* d = c.decay({54, 135, 0});
  ASSERT_NE(d, nullptr);
  EXPECT_NEAR(d->decayConstant, kLn2 / 100.0, 1e-15);
}

TEST(Chain, StableNuclideHasZeroDecayConstant) {
  DepletionChain c;
  c.setDecay({2, 4, 0},
             DecayData{.halfLife = 0.0, .decayConstant = 999.0});  // bogus input constant
  EXPECT_EQ(c.decay({2, 4, 0})->decayConstant, 0.0);
  c.setDecay({1, 1, 0}, DecayData{.halfLife = std::numeric_limits<double>::infinity()});
  EXPECT_EQ(c.decay({1, 1, 0})->decayConstant, 0.0);
}

TEST(Chain, DecayLookupAbsentReturnsNull) {
  DepletionChain c;
  EXPECT_EQ(c.decay({92, 238, 0}), nullptr);
}

// ---------------------------------------------------------------------------
// Fission yields
// ---------------------------------------------------------------------------
TEST(Chain, FissionYieldsRegisterParentAndProducts) {
  DepletionChain c;
  FissionYields y{.energy = 0.0253, .products = {{{54, 135, 0}, 0.064}, {{42, 99, 0}, 0.061}}};
  c.addFissionYields({92, 235, 0}, y);
  EXPECT_GE(c.indexOf({92, 235, 0}), 0);
  EXPECT_GE(c.indexOf({54, 135, 0}), 0);
  EXPECT_GE(c.indexOf({42, 99, 0}), 0);
}

TEST(Chain, NearestYieldsByEnergy) {
  DepletionChain c;
  FissionYields thermal{.energy = 0.0253, .products = {{{54, 135, 0}, 0.06}}};
  FissionYields fast{.energy = 5.0e5, .products = {{{54, 135, 0}, 0.05}}};
  c.addFissionYields({92, 235, 0}, thermal);
  c.addFissionYields({92, 235, 0}, fast);

  EXPECT_DOUBLE_EQ(c.nearestYields({92, 235, 0}, 1.0)->energy, 0.0253);
  EXPECT_DOUBLE_EQ(c.nearestYields({92, 235, 0}, 1.0e6)->energy, 5.0e5);
  EXPECT_EQ(c.nearestYields({94, 239, 0}, 1.0), nullptr);
}

// ---------------------------------------------------------------------------
// Decay matrix assembly
// ---------------------------------------------------------------------------
TEST(DecayMatrix, SingleDecayEntries) {
  DepletionChain c;
  const Zai A{53, 135, 0}, B{54, 135, 0};
  c.add(A);
  c.add(B);
  c.setDecay(
      A, DecayData{
             .halfLife = 100.0,
             .modes = {DecayMode{
                 .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});  // A -> B

  auto M = c.decayMatrix();
  const double lam = kLn2 / 100.0;
  int ia = c.indexOf(A), ib = c.indexOf(B);
  EXPECT_NEAR(M.coeff(ia, ia), -lam, 1e-15);  // removal
  EXPECT_NEAR(M.coeff(ib, ia), +lam, 1e-15);  // production of daughter
  EXPECT_NEAR(M.coeff(ib, ib), 0.0, 1e-15);   // B stable
}

TEST(DecayMatrix, BranchingSplitsByRatio) {
  DepletionChain c;
  const Zai A{1, 100, 0}, B{2, 100, 0}, D{1, 99, 0};
  c.add(A);
  c.add(B);
  c.add(D);
  // beta- (Z+1 -> B) with BR 0.7 ; neutron emission (A-1 -> D) with BR 0.3
  c.setDecay(
      A, DecayData{
             .halfLife = 10.0,
             .modes = {
                 DecayMode{.rtyp = 1.0, .branching = 0.7, .finalState = 0, .isFission = false},
                 DecayMode{.rtyp = 5.0, .branching = 0.3, .finalState = 0, .isFission = false}}});

  auto M = c.decayMatrix();
  const double lam = kLn2 / 10.0;
  int ia = c.indexOf(A);
  EXPECT_NEAR(M.coeff(ia, ia), -lam, 1e-15);
  EXPECT_NEAR(M.coeff(c.indexOf(B), ia), 0.7 * lam, 1e-15);
  EXPECT_NEAR(M.coeff(c.indexOf(D), ia), 0.3 * lam, 1e-15);
}

TEST(DecayMatrix, ConservativeChainHasZeroColumnSums) {
  // Every tracked decay maps one nucleus to exactly one tracked nucleus, so
  // each column must sum to zero (atom-number conserving).
  DepletionChain c;
  const Zai A{53, 135, 0}, B{54, 135, 0}, D{55, 135, 0};
  c.add(A);
  c.add(B);
  c.add(D);
  c.setDecay(A,
             DecayData{.halfLife = 100.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  c.setDecay(B,
             DecayData{.halfLife = 300.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});

  auto M = c.decayMatrix();
  for (int col = 0; col < M.cols(); ++col) {
    double s = 0.0;
    for (int row = 0; row < M.rows(); ++row)
      s += M.coeff(row, col);
    EXPECT_NEAR(s, 0.0, 1e-14) << "column " << col;
  }
}

// close() registers every reachable daughter so production is never dropped. Here
// Cs-135 (Xe-135's daughter) is absent until close() adds it, after which the decay
// matrix conserves atoms (zero column sums). Idempotent on an already-closed chain.
TEST(Chain, CloseRegistersDaughtersAndRestoresConservation) {
  DepletionChain c;
  const Zai A{53, 135, 0}, B{54, 135, 0};  // I-135 -> Xe-135 -> Cs-135
  c.setDecay(A,
             DecayData{.halfLife = 100.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  c.setDecay(B,
             DecayData{.halfLife = 300.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  EXPECT_EQ(c.indexOf({55, 135, 0}), -1);  // daughter absent -> would be dropped

  EXPECT_EQ(c.close(), 1);  // Cs-135 registered
  EXPECT_GE(c.indexOf({55, 135, 0}), 0);
  EXPECT_EQ(c.close(), 0);  // idempotent

  auto M = c.decayMatrix();
  for (int col = 0; col < M.cols(); ++col) {
    double s = 0.0;
    for (int row = 0; row < M.rows(); ++row)
      s += M.coeff(row, col);
    EXPECT_NEAR(s, 0.0, 1e-14) << "column " << col;
  }
}

TEST(DecayMatrix, SpontaneousFissionFeedsYields) {
  DepletionChain c;
  const Zai Cf{98, 252, 0}, P1{54, 140, 0}, P2{44, 108, 0};
  c.add(Cf);
  FissionYields sfy{.energy = 0.0, .products = {{P1, 1.2}, {P2, 0.8}}};
  c.addFissionYields(Cf, sfy);
  // single SF branch, branching 1.0
  c.setDecay(Cf,
             DecayData{.halfLife = 8.3e8,
                       .modes = {DecayMode{
                           .rtyp = 6.0, .branching = 1.0, .finalState = 0, .isFission = true}}});

  auto M = c.decayMatrix();
  const double lam = kLn2 / 8.3e8;
  int icf = c.indexOf(Cf);
  EXPECT_NEAR(M.coeff(icf, icf), -lam, 1e-20);
  EXPECT_NEAR(M.coeff(c.indexOf(P1), icf), 1.2 * lam, 1e-20);
  EXPECT_NEAR(M.coeff(c.indexOf(P2), icf), 0.8 * lam, 1e-20);
}

// ---------------------------------------------------------------------------
// Reaction / fission source helpers
// ---------------------------------------------------------------------------
TEST(Reactions, AddReactionRemovesParentProducesProduct) {
  DepletionChain c;
  const Zai A{92, 235, 0}, B{92, 236, 0};
  c.add(A);
  c.add(B);
  std::vector<Eigen::Triplet<double>> t;
  c.addReaction(t, A, B, 3.0);  // (n,gamma)
  auto M = c.finalize(t);
  EXPECT_NEAR(M.coeff(c.indexOf(A), c.indexOf(A)), -3.0, 1e-15);
  EXPECT_NEAR(M.coeff(c.indexOf(B), c.indexOf(A)), +3.0, 1e-15);
}

TEST(Reactions, FissionSourceWeightsByYield) {
  DepletionChain c;
  const Zai U{92, 235, 0}, P1{54, 135, 0}, P2{42, 99, 0};
  c.add(U);
  FissionYields y{.energy = 0.0253, .products = {{P1, 0.064}, {P2, 0.061}}};
  c.addFissionYields(U, y);

  std::vector<Eigen::Triplet<double>> t;
  c.addFissionSource(t, U, /*fissionRate=*/2.0, /*energy=*/0.0253);
  auto M = c.finalize(t);
  EXPECT_NEAR(M.coeff(c.indexOf(U), c.indexOf(U)), -2.0, 1e-15);
  EXPECT_NEAR(M.coeff(c.indexOf(P1), c.indexOf(U)), 2.0 * 0.064, 1e-15);
  EXPECT_NEAR(M.coeff(c.indexOf(P2), c.indexOf(U)), 2.0 * 0.061, 1e-15);
}

TEST(Reactions, AddReactionEarlyReturns) {
  DepletionChain c;
  const Zai A{92, 235, 0}, B{92, 236, 0}, X{99, 999, 0};  // X never registered
  c.add(A);
  c.add(B);
  std::vector<Eigen::Triplet<double>> t;

  c.addReaction(t, A, X, 1.0);  // unregistered product
  c.addReaction(t, X, A, 1.0);  // unregistered parent
  c.addReaction(t, A, B, 0.0);  // zero rate
  EXPECT_TRUE(t.empty());
}

TEST(Reactions, AddFissionSourceEarlyReturns) {
  DepletionChain c;
  const Zai U{92, 235, 0}, Z{99, 999, 0};
  c.add(U);
  std::vector<Eigen::Triplet<double>> t;

  c.addFissionSource(t, Z, 1.0, 0.0253);  // unregistered parent
  c.addFissionSource(t, U, 0.0, 0.0253);  // zero rate
  c.addFissionSource(t, U, 1.0, 0.0253);  // registered, but no yields supplied
  EXPECT_TRUE(t.empty());
}

TEST(Reactions, FinalizeSumsDuplicateTriplets) {
  DepletionChain c;
  c.add({1, 1, 0});
  std::vector<Eigen::Triplet<double>> t{{0, 0, 1.5}, {0, 0, 2.5}};
  auto M = c.finalize(t);
  EXPECT_NEAR(M.coeff(0, 0), 4.0, 1e-15);
}

// When built without -DWITH_ENDFTK, the readers are no-ops that load nothing.
#ifndef WITH_ENDFTK
TEST(EndfReader, StubsLoadNothingWithoutEndftk) {
  DepletionChain c;
  EXPECT_EQ(loadDecayData(c, "ignored.endf"), 0);
  EXPECT_EQ(loadFissionYields(c, "ignored.endf"), 0);
  EXPECT_EQ(c.size(), 0);
}
#endif

// Spontaneous fission with no SFY table registered: the parent is still
// removed, but no products are produced (column sum is negative).
TEST(DecayMatrix, SpontaneousFissionWithoutYieldsJustRemoves) {
  DepletionChain c;
  const Zai Cm{96, 244, 0};
  c.add(Cm);
  c.setDecay(Cm,
             DecayData{.halfLife = 5.7e8,
                       .modes = {DecayMode{
                           .rtyp = 6.0, .branching = 1.0, .finalState = 0, .isFission = true}}});
  auto M = c.decayMatrix();
  const double lam = kLn2 / 5.7e8;
  int i = c.indexOf(Cm);
  EXPECT_NEAR(M.coeff(i, i), -lam, 1e-20);
  double col = 0.0;
  for (int r = 0; r < M.rows(); ++r)
    col += M.coeff(r, i);
  EXPECT_NEAR(col, -lam, 1e-20);  // mass leaves the tracked system
}

// An unregistered daughter means the matrix silently loses that production and
// stops conserving atoms. Callers that want to check for it should not have to
// scrape stderr.
TEST(DecayMatrix, ReportsDroppedDaughtersThroughTheOutParameter) {
  DepletionChain c;
  c.setDecay({40, 97, 0},
             DecayData{.halfLife = 1000.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  c.setDecay({52, 133, 0},
             DecayData{.halfLife = 900.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});

  int dropped = -1;
  auto before = c.decayMatrix(&dropped);
  EXPECT_EQ(dropped, 2) << "neither daughter is registered yet";
  EXPECT_EQ(before.nonZeros(), 2) << "only the two removal terms survive";

  EXPECT_EQ(c.close(), 2);
  dropped = -1;
  auto after = c.decayMatrix(&dropped);
  EXPECT_EQ(dropped, 0) << "close() registers every reachable daughter";

  // With the chain closed, every column sums to zero and atoms are conserved.
  for (int col = 0; col < after.cols(); ++col) {
    double sum = 0.0;
    for (int row = 0; row < after.rows(); ++row)
      sum += after.coeff(row, col);
    EXPECT_NEAR(sum, 0.0, 1e-18) << "column " << col;
  }
}

// A nuclide that was registered with setDecay() but is stable (halfLife=0)
// must be skipped during matrix assembly without producing a diagonal entry.
TEST(DecayMatrix, ExplicitlyStableNuclideContributesNothing) {
  DepletionChain c;
  const Zai S{2, 4, 0};  // 4He, stable, registered via setDecay
  c.setDecay(S, DecayData{.halfLife = 0.0});
  auto M = c.decayMatrix();
  EXPECT_EQ(M.nonZeros(), 0);
}

// A non-fission mode whose RTYP sequence still encodes fission (e.g. beta-
// followed by spontaneous fission, RTYP 1.6) must drop the production term
// rather than try to apply a meaningless daughter.
TEST(DecayMatrix, FissionInMultiStepRtypIsSkipped) {
  DepletionChain c;
  const Zai parent{92, 238, 0}, dummy{93, 238, 0};
  c.add(parent);
  c.add(dummy);
  // isFission=false but RTYP contains code 6 -> applyDecay flags fission
  // and the production triplet must be omitted.
  c.setDecay(parent,
             DecayData{.halfLife = 1.0e9,
                       .modes = {DecayMode{
                           .rtyp = 1.6, .branching = 1.0, .finalState = 0, .isFission = false}}});
  auto M = c.decayMatrix();
  EXPECT_NEAR(M.coeff(c.indexOf(parent), c.indexOf(parent)), -(kLn2 / 1.0e9), 1e-20);
  // No production of the would-be daughter, since fission was encoded.
  EXPECT_EQ(M.nonZeros(), 1);
}

// A daughter that was never registered in the chain has its production term
// silently dropped (documented behaviour: register reachable daughters).
TEST(DecayMatrix, UnregisteredDaughterIsDropped) {
  DepletionChain c;
  const Zai A{53, 135, 0};  // beta- -> Xe-135, which we deliberately omit
  c.add(A);
  c.setDecay(A,
             DecayData{.halfLife = 100.0,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  auto M = c.decayMatrix();
  EXPECT_EQ(c.indexOf({54, 135, 0}), -1);              // daughter absent
  EXPECT_NEAR(M.coeff(0, 0), -(kLn2 / 100.0), 1e-15);  // parent still removed
  EXPECT_EQ(M.nonZeros(), 1);                          // no production term added
}

// ---------------------------------------------------------------------------
// Explicit decay daughters
// ---------------------------------------------------------------------------
// A source that names the product of a decay mode (an OpenMC depletion chain
// does) must have that product honored even where the RTYP transition rules
// would derive a different one.
TEST(DecayMatrix, ExplicitDaughterOverridesRtypDerivedOne) {
  DepletionChain c;
  const Zai A{1, 100, 0}, derived{2, 100, 0}, named{3, 100, 0};
  c.add(A);
  c.add(derived);
  c.add(named);
  // rtyp 1.0 (beta-) would produce `derived`; the explicit daughter says `named`.
  c.setDecay(A, DecayData{.halfLife = 10.0,
                          .modes = {DecayMode{.rtyp = 1.0,
                                              .branching = 1.0,
                                              .finalState = 0,
                                              .isFission = false,
                                              .daughter = named}}});
  auto M = c.decayMatrix();
  const double lam = kLn2 / 10.0;
  EXPECT_NEAR(M.coeff(c.indexOf(A), c.indexOf(A)), -lam, 1e-15);
  EXPECT_NEAR(M.coeff(c.indexOf(named), c.indexOf(A)), lam, 1e-15);
  EXPECT_EQ(M.coeff(c.indexOf(derived), c.indexOf(A)), 0.0);
}

// close() must register the same daughter decayTriplets() will produce, so an
// explicit daughter that is not yet in the chain is what close() adds.
TEST(Chain, CloseRegistersExplicitDaughter) {
  DepletionChain c;
  const Zai A{1, 100, 0}, named{3, 100, 0};
  c.setDecay(A, DecayData{.halfLife = 10.0,
                          .modes = {DecayMode{.rtyp = 1.0,
                                              .branching = 1.0,
                                              .finalState = 0,
                                              .isFission = false,
                                              .daughter = named}}});
  EXPECT_EQ(c.indexOf(named), -1);
  EXPECT_EQ(c.indexOf({2, 100, 0}), -1);

  EXPECT_EQ(c.close(), 1);
  EXPECT_GE(c.indexOf(named), 0);
  EXPECT_EQ(c.indexOf({2, 100, 0}), -1) << "the rtyp-derived daughter is not registered";

  int dropped = -1;
  auto M = c.decayMatrix(&dropped);
  EXPECT_EQ(dropped, 0);
  for (int col = 0; col < M.cols(); ++col) {
    double sum = 0.0;
    for (int row = 0; row < M.rows(); ++row)
      sum += M.coeff(row, col);
    EXPECT_NEAR(sum, 0.0, 1e-18) << "column " << col;
  }
}

// An explicit daughter on a fission mode is meaningless and must be ignored:
// the products come from the yield table, and close() must not register it.
TEST(DecayMatrix, ExplicitDaughterIgnoredForFissionMode) {
  DepletionChain c;
  const Zai Cm{96, 244, 0}, bogus{1, 1, 0};
  c.add(Cm);
  c.setDecay(Cm, DecayData{.halfLife = 5.7e8,
                           .modes = {DecayMode{.rtyp = 6.0,
                                               .branching = 1.0,
                                               .finalState = 0,
                                               .isFission = true,
                                               .daughter = bogus}}});
  EXPECT_EQ(c.close(), 0);
  EXPECT_EQ(c.indexOf(bogus), -1);
  auto M = c.decayMatrix();
  EXPECT_EQ(M.nonZeros(), 1);  // removal only; no SFY table supplied
}
