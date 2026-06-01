#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include "cram/chain.hpp"
#include "cram/chain_xml.hpp"
#include "cram/deplete.hpp"

using namespace cram;

namespace {

// A trimmed depletion_chain in the OpenMC schema, covering decay (alpha + sf),
// the reaction types we model, a "Nothing" target, and a fission-yield block.
const char* kChainXml = R"XML(<?xml version="1.0"?>
<depletion_chain>
  <nuclide name="U235" half_life="2.22102e16" decay_modes="2" reactions="3">
    <decay type="sf" target="U235" branching_ratio="7.0e-11"/>
    <decay type="alpha" target="Th231" branching_ratio="0.99999999993"/>
    <reaction type="(n,gamma)" Q="6545200.0" target="U236"/>
    <reaction type="(n,2n)" Q="-5297781.0" target="U234"/>
    <reaction type="(n,3n)" Q="-12142300.0" target="U233"/>
    <reaction type="(n,4n)" Q="-17885600.0" target="U232"/>
    <reaction type="(n,a)" Q="-1000000.0" target="Th232"/>
    <reaction type="(n,d)" Q="-2000000.0" target="Nothing"/>
    <reaction type="fission" Q="193405400.0"/>
    <neutron_fission_yields>
      <energies>0.0253</energies>
      <fission_yields energy="0.0253">
        <products>Xe135 I135 Zr95</products>
        <data>0.00254 0.0292 0.0653</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="U238" half_life="1.409967e17" decay_modes="1" reactions="2">
    <decay type="alpha" target="Th234" branching_ratio="1.0"/>
    <reaction type="(n,gamma)" Q="4806800.0" target="U239"/>
    <reaction type="(n,p)" Q="-1234000.0" target="Nothing"/>
  </nuclide>
  <nuclide name="Th231" half_life="91872.0" decay_modes="1" reactions="0">
    <decay type="beta-" target="Pa231" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="Xe135" half_life="32904.0" decay_modes="1" reactions="0">
    <decay type="beta-" target="Cs135" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="I135" half_life="23652.0" decay_modes="1" reactions="0">
    <decay type="beta-" target="Xe135" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="Gd157" reactions="1">
    <reaction type="(n,gamma)" Q="7937400.0" target="Gd158"/>
  </nuclide>
  <nuclide name="Ag109" reactions="2">
    <reaction type="(n,gamma)" Q="6809000.0" target="Ag110" branching_ratio="0.954"/>
    <reaction type="(n,gamma)" Q="6809000.0" target="Ag110_m1" branching_ratio="0.046"/>
  </nuclide>
</depletion_chain>
)XML";

std::string writeTempChain() {
  std::string path = std::string(std::tmpnam(nullptr)) + "_chain.xml";
  std::ofstream(path) << kChainXml;
  return path;
}

}  // namespace

TEST(ChainXml, ParseNuclideName) {
  Zai z;
  EXPECT_TRUE(parseNuclideName("U235", z));
  EXPECT_EQ(z, (Zai{92, 235, 0}));
  EXPECT_TRUE(parseNuclideName("Am242_m1", z));
  EXPECT_EQ(z, (Zai{95, 242, 1}));
  EXPECT_TRUE(parseNuclideName("Th231", z));
  EXPECT_EQ(z, (Zai{90, 231, 0}));
  EXPECT_TRUE(parseNuclideName("Xe135", z));
  EXPECT_EQ(z, (Zai{54, 135, 0}));

  EXPECT_FALSE(parseNuclideName("Zz999", z));  // bad element
  EXPECT_FALSE(parseNuclideName("U", z));      // no mass
  EXPECT_FALSE(parseNuclideName("235", z));    // no symbol
}

TEST(ChainXml, ElementSymbol) {
  EXPECT_EQ(elementSymbol(92), "U");
  EXPECT_EQ(elementSymbol(54), "Xe");
  EXPECT_EQ(elementSymbol(1), "H");
  EXPECT_EQ(elementSymbol(0), "");
  EXPECT_EQ(elementSymbol(200), "");
}

TEST(ChainXml, LoadChain) {
  const std::string path = writeTempChain();
  DepletionChain chain;
  auto reactions = loadDepletionChainXml(chain, path);
  std::remove(path.c_str());

  // Nuclides registered (including decay/reaction targets).
  EXPECT_GE(chain.indexOf({92, 235, 0}), 0);  // U235
  EXPECT_GE(chain.indexOf({90, 231, 0}), 0);  // Th231 (alpha daughter)
  EXPECT_GE(chain.indexOf({54, 135, 0}), 0);  // Xe135

  // Decay data parsed: half-life and an alpha branch to Th231.
  const DecayData* u235 = chain.decay({92, 235, 0});
  ASSERT_NE(u235, nullptr);
  EXPECT_NEAR(u235->halfLife, 2.22102e16, 1e9);
  bool hasAlphaToTh231 = false, hasSf = false;
  for (const auto& m : u235->modes) {
    if (m.isFission)
      hasSf = true;
    if (m.hasDaughter && m.daughter == Zai{90, 231, 0})
      hasAlphaToTh231 = true;
  }
  EXPECT_TRUE(hasAlphaToTh231);
  EXPECT_TRUE(hasSf);

  // Fission yields parsed for U235.
  const FissionYields* y = chain.nearestYields({92, 235, 0}, 0.0253);
  ASSERT_NE(y, nullptr);
  EXPECT_EQ(y->products.size(), 3u);

  // Reaction topology: U235 fission + (n,gamma)->U236 + (n,2n)->U234; U238
  // (n,gamma)->U239 + (n,p)->Nothing (no tracked target).
  int nFission = 0, nGammaToU236 = 0, nNothing = 0;
  for (const auto& r : reactions) {
    if (r.parent == Zai{92, 235, 0} && r.type == ReactionType::Fission)
      ++nFission;
    if (r.parent == Zai{92, 235, 0} && r.type == ReactionType::NGamma && r.hasTarget &&
        r.target == Zai{92, 236, 0})
      ++nGammaToU236;
    if (r.parent == Zai{92, 238, 0} && r.type == ReactionType::NProton && !r.hasTarget)
      ++nNothing;
  }
  EXPECT_EQ(nFission, 1);
  EXPECT_EQ(nGammaToU236, 1);
  EXPECT_EQ(nNothing, 1);
}

TEST(ChainXml, StableNuclideHasNoDecayAndUnknownReactionsAreSkipped) {
  const std::string path = writeTempChain();
  DepletionChain chain;
  auto reactions = loadDepletionChainXml(chain, path);
  std::remove(path.c_str());

  // Gd157 has no half_life -> registered but no decay data.
  EXPECT_GE(chain.indexOf({64, 157, 0}), 0);
  EXPECT_EQ(chain.decay({64, 157, 0}), nullptr);

  // The unknown "(n,d)" reaction type is dropped; the modeled types survive.
  for (const auto& r : reactions)
    EXPECT_TRUE(r.type == ReactionType::Fission || r.type == ReactionType::NGamma ||
                r.type == ReactionType::N2n || r.type == ReactionType::N3n ||
                r.type == ReactionType::N4n || r.type == ReactionType::NAlpha ||
                r.type == ReactionType::NProton);
  bool sawN3n = false, sawN4n = false, sawAlpha = false;
  for (const auto& r : reactions) {
    sawN3n |= (r.type == ReactionType::N3n);
    sawN4n |= (r.type == ReactionType::N4n);
    sawAlpha |= (r.type == ReactionType::NAlpha);
  }
  EXPECT_TRUE(sawN3n && sawN4n && sawAlpha);

  // Branched (n,gamma): Ag109 -> Ag110 (0.954) + Ag110_m1 (0.046).
  double brGround = 0.0, brMeta = 0.0;
  for (const auto& r : reactions) {
    if (r.parent != Zai{47, 109, 0} || r.type != ReactionType::NGamma)
      continue;
    if (r.target == Zai{47, 110, 0})
      brGround = r.branching;
    else if (r.target == Zai{47, 110, 1})
      brMeta = r.branching;
  }
  EXPECT_DOUBLE_EQ(brGround, 0.954);
  EXPECT_DOUBLE_EQ(brMeta, 0.046);
}

TEST(ChainXml, MalformedAndMissingRootThrow) {
  DepletionChain chain;
  EXPECT_THROW(loadDepletionChainXml(chain, "/no/such/file_xyz.xml"), std::runtime_error);

  std::string bad = std::string(std::tmpnam(nullptr)) + "_bad.xml";
  std::ofstream(bad) << "<not_a_chain><foo/></not_a_chain>";
  EXPECT_THROW(loadDepletionChainXml(chain, bad), std::runtime_error);
  std::remove(bad.c_str());
}

// A reaction whose product is not tracked (here U238 (n,p) -> "Nothing") must
// still remove the parent under irradiation.
TEST(ChainXml, UntrackedReactionProductStillRemovesParent) {
  const std::string path = writeTempChain();
  DepletionChain chain;
  auto reactions = loadDepletionChainXml(chain, path);
  std::remove(path.c_str());

  DepletionSystem sys(chain);
  // Give U238 only its (n,p)->Nothing channel a cross section.
  sys.setReactions({92, 238, 0}, {{ReactionType::NProton, {}, 5.0}});
  sys.setConstantFlux(1.0e16);

  Eigen::SparseMatrix<double> A = sys.assemble(Eigen::VectorXd::Zero(chain.size()));
  const int i = chain.indexOf({92, 238, 0});
  EXPECT_LT(A.coeff(i, i), 0.0);  // parent removed even though product untracked
}
