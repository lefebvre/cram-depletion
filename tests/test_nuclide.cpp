#include <gtest/gtest.h>

#include <cmath>

#include "cram/nuclide.hpp"

using namespace cram;

// ---------------------------------------------------------------------------
// Zai identity
// ---------------------------------------------------------------------------
TEST(Zai, KeyPacking) {
  EXPECT_EQ((Zai{92, 235, 0}).key(), 922350);
  EXPECT_EQ((Zai{92, 235, 1}).key(), 922351);
  EXPECT_EQ((Zai{1, 1, 0}).key(), 10010);
  EXPECT_EQ((Zai{0, 1, 0}).key(), 10);  // a neutron pseudo-nuclide
}

TEST(Zai, Equality) {
  EXPECT_EQ((Zai{54, 135, 0}), (Zai{54, 135, 0}));
  EXPECT_FALSE((Zai{54, 135, 0}) == (Zai{54, 135, 1}));
  EXPECT_FALSE((Zai{54, 135, 0}) == (Zai{55, 135, 0}));
}

TEST(Zai, StringForm) {
  EXPECT_EQ((Zai{53, 135, 0}).str(), "53-135");
  EXPECT_EQ((Zai{43, 99, 1}).str(), "43-99m1");
}

TEST(Zai, HashUsableInMap) {
  std::unordered_map<Zai, int, ZaiHash> m;
  m[{92, 238, 0}] = 7;
  m[{92, 238, 1}] = 8;  // distinct metastable state -> distinct key
  EXPECT_EQ(m.size(), 2u);
  EXPECT_EQ((m[{92, 238, 0}]), 7);
}

// ---------------------------------------------------------------------------
// RTYP decode
// ---------------------------------------------------------------------------
namespace {
// decayParticleSequence returns a fixed-capacity sequence rather than a vector;
// materialise it so the expectations below stay readable.
std::vector<int> codes(double rtyp) {
  const auto seq = decayParticleSequence(rtyp);
  return std::vector<int>(seq.begin(), seq.end());
}
}  // namespace

TEST(DecaySequence, SingleParticle) {
  EXPECT_EQ(codes(1.0), (std::vector<int>{1}));
  EXPECT_EQ(codes(4.0), (std::vector<int>{4}));
  EXPECT_EQ(codes(6.0), (std::vector<int>{6}));
}

TEST(DecaySequence, MultiParticle) {
  EXPECT_EQ(codes(1.5), (std::vector<int>{1, 5}));  // beta- , delayed n
  EXPECT_EQ(codes(1.4), (std::vector<int>{1, 4}));  // beta- , alpha
  EXPECT_EQ(codes(2.7), (std::vector<int>{2, 7}));  // EC , proton
  EXPECT_EQ(codes(1.55), (std::vector<int>{1, 5, 5}));
}

TEST(DecaySequence, NeverExceedsItsFixedCapacity) {
  // 1.23456789 decodes to the integer part plus all eight fractional digits,
  // which is the longest sequence RTYP can encode at kScale = 1e8.
  EXPECT_EQ(decayParticleSequence(1.23456789).size(), 9);
}

TEST(DecaySequence, NegativeRtypIsEmpty) {
  EXPECT_TRUE(decayParticleSequence(-1.0).empty());
  EXPECT_TRUE(decayParticleSequence(std::nan("")).empty());
}

// ---------------------------------------------------------------------------
// applyDecay: daughter computation
// ---------------------------------------------------------------------------
namespace {
Zai daughter(Zai parent, double rtyp, int rfs = 0) {
  bool fission = false;
  return applyDecay(parent, rtyp, rfs, fission);
}
bool isFission(Zai parent, double rtyp) {
  bool fission = false;
  applyDecay(parent, rtyp, 0, fission);
  return fission;
}
}  // namespace

TEST(ApplyDecay, BetaMinus) {  // n -> p, Z+1
  EXPECT_EQ(daughter({6, 14, 0}, 1.0), (Zai{7, 14, 0}));
  EXPECT_EQ(daughter({53, 135, 0}, 1.0), (Zai{54, 135, 0}));
}

TEST(ApplyDecay, ElectronCapture) {  // p -> n, Z-1
  EXPECT_EQ(daughter({11, 22, 0}, 2.0), (Zai{10, 22, 0}));
}

TEST(ApplyDecay, Alpha) {  // Z-2, A-4
  EXPECT_EQ(daughter({92, 238, 0}, 4.0), (Zai{90, 234, 0}));
  EXPECT_EQ(daughter({94, 239, 0}, 4.0), (Zai{92, 235, 0}));
}

TEST(ApplyDecay, IsomericTransition) {  // no Z/A change; lands on requested state
  EXPECT_EQ(daughter({43, 99, 1}, 3.0, 0), (Zai{43, 99, 0}));
}

TEST(ApplyDecay, NeutronAndProtonEmission) {
  EXPECT_EQ(daughter({8, 17, 0}, 5.0), (Zai{8, 16, 0}));  // neutron emission
  EXPECT_EQ(daughter({7, 12, 0}, 7.0), (Zai{6, 11, 0}));  // proton emission
}

TEST(ApplyDecay, MultiStepChains) {
  EXPECT_EQ(daughter({53, 137, 0}, 1.5), (Zai{54, 136, 0}));  // beta- then delayed n
  EXPECT_EQ(daughter({3, 8, 0}, 1.4), (Zai{2, 4, 0}));        // beta- then alpha
}

TEST(ApplyDecay, FinalStateApplied) {
  EXPECT_EQ(daughter({27, 60, 0}, 1.0, 1), (Zai{28, 60, 1}));  // beta- to a metastable daughter
}

TEST(ApplyDecay, SpontaneousFissionFlagged) {
  EXPECT_TRUE(isFission({98, 252, 0}, 6.0));
  EXPECT_FALSE(isFission({92, 238, 0}, 4.0));
}

TEST(ApplyDecay, UnknownCodeLeavesNuclideUnchanged) {
  // ENDF code 10 ("unknown"): no Z/A change, but the requested final
  // isomeric state is still applied.
  bool fission = false;
  Zai d = applyDecay({50, 120, 0}, 10.0, 1, fission);
  EXPECT_FALSE(fission);
  EXPECT_EQ(d, (Zai{50, 120, 1}));
}
