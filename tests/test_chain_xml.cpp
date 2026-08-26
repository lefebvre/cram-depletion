#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
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

// Entries the reader cannot use: a nuclide name, a decay target and a yield
// product that do not parse, and a reaction type outside the modeled set.
const char* kDirtyChainXml = R"XML(<depletion_chain>
  <nuclide name="Zz999" half_life="1.0"/>
  <nuclide name="I135" half_life="23652.0" decay_modes="2" reactions="1">
    <decay type="beta-" target="Xe135" branching_ratio="0.9"/>
    <decay type="beta-" target="Qq135" branching_ratio="0.1"/>
    <reaction type="(n,t)" Q="0.0" target="Nothing"/>
    <neutron_fission_yields>
      <fission_yields energy="0.0253">
        <products>Xe135 Bogus1</products>
        <data>0.5 0.5</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="Xe135"/>
</depletion_chain>)XML";

std::optional<DecayMode> modeTo(const DecayData& d, const Zai& daughter) {
  for (const auto& m : d.modes)
    if (m.daughter == daughter)
      return m;
  return std::nullopt;
}

// A temporary file that is removed when the test scope ends.
struct TempFile {
  std::filesystem::path path;
  explicit TempFile(const char* stem, const char* text)
      : path(std::filesystem::temp_directory_path() / stem) {
    std::ofstream(path) << text;
  }
  ~TempFile() { std::filesystem::remove(path); }
};

}  // namespace

TEST(ChainXml, ParseNuclideName) {
  EXPECT_EQ(parseNuclideName("U235"), (Zai{92, 235, 0}));
  EXPECT_EQ(parseNuclideName("Am242_m1"), (Zai{95, 242, 1}));
  EXPECT_EQ(parseNuclideName("Th231"), (Zai{90, 231, 0}));
  EXPECT_EQ(parseNuclideName("Xe135"), (Zai{54, 135, 0}));

  EXPECT_FALSE(parseNuclideName("Zz999").has_value());  // bad element
  EXPECT_FALSE(parseNuclideName("U").has_value());      // no mass
  EXPECT_FALSE(parseNuclideName("235").has_value());    // no symbol
  EXPECT_FALSE(parseNuclideName("").has_value());
}

TEST(ChainXml, ElementSymbol) {
  EXPECT_EQ(elementSymbol(92), "U");
  EXPECT_EQ(elementSymbol(54), "Xe");
  EXPECT_EQ(elementSymbol(1), "H");
  EXPECT_EQ(elementSymbol(0), "");
  EXPECT_EQ(elementSymbol(200), "");
}

TEST(ChainXml, LoadChain) {
  DepletionChain chain;
  ChainXmlDiagnostics diag;
  auto reactions = loadDepletionChainXmlString(chain, kChainXml, &diag);
  EXPECT_EQ(diag.unparsedNuclides, 0);
  EXPECT_EQ(diag.unparsedDecayTargets, 0);
  EXPECT_EQ(diag.unparsedYieldProducts, 0);
  EXPECT_EQ(diag.unmodeledReactions, 1);  // the (n,d) entry

  // Nuclides registered (including decay/reaction targets).
  EXPECT_GE(chain.indexOf({92, 235, 0}), 0);  // U235
  EXPECT_GE(chain.indexOf({90, 231, 0}), 0);  // Th231 (alpha daughter)
  EXPECT_GE(chain.indexOf({54, 135, 0}), 0);  // Xe135

  // Decay data parsed: half-life, an alpha branch to Th231 carried as an
  // explicit daughter, and an SF branch.
  const DecayData* u235 = chain.decay({92, 235, 0});
  ASSERT_NE(u235, nullptr);
  EXPECT_NEAR(u235->halfLife, 2.22102e16, 1e9);
  const auto alpha = modeTo(*u235, Zai{90, 231, 0});
  ASSERT_TRUE(alpha.has_value());
  EXPECT_DOUBLE_EQ(alpha->branching, 0.99999999993);
  EXPECT_FALSE(alpha->isFission);
  bool hasSf = false;
  for (const auto& m : u235->modes)
    hasSf |= m.isFission;
  EXPECT_TRUE(hasSf);

  // Fission yields parsed for U235.
  const FissionYields* y = chain.nearestYields({92, 235, 0}, 0.0253);
  ASSERT_NE(y, nullptr);
  EXPECT_EQ(y->products.size(), 3u);

  // Reaction topology: U235 fission + (n,gamma)->U236 + (n,2n)->U234; U238
  // (n,gamma)->U239 + (n,p)->Nothing (no tracked target).
  int nFission = 0, nGammaToU236 = 0, nNothing = 0;
  for (const auto& r : reactions) {
    if (r.parent == Zai{92, 235, 0} && r.type == ReactionType::Fission) {
      ++nFission;
      EXPECT_DOUBLE_EQ(r.q, 193405400.0);
    }
    if (r.parent == Zai{92, 235, 0} && r.type == ReactionType::NGamma &&
        r.target == Zai{92, 236, 0})
      ++nGammaToU236;
    if (r.parent == Zai{92, 238, 0} && r.type == ReactionType::NProton && !r.target)
      ++nNothing;
  }
  EXPECT_EQ(nFission, 1);
  EXPECT_EQ(nGammaToU236, 1);
  EXPECT_EQ(nNothing, 1);
}

// The explicit targets route the decay matrix: Th231 -> Pa231 as the file says,
// and the chain closes without dropping any daughter.
TEST(ChainXml, DecayMatrixFollowsExplicitTargets) {
  DepletionChain chain;
  loadDepletionChainXmlString(chain, kChainXml);
  chain.close();
  int dropped = -1;
  auto M = chain.decayMatrix(&dropped);
  EXPECT_EQ(dropped, 0);
  const int th = chain.indexOf({90, 231, 0});
  const int pa = chain.indexOf({91, 231, 0});
  ASSERT_GE(pa, 0);
  EXPECT_GT(M.coeff(pa, th), 0.0);
  EXPECT_NEAR(M.coeff(pa, th), -M.coeff(th, th), 1e-25);
}

TEST(ChainXml, StableNuclideHasNoDecayAndUnknownReactionsAreSkipped) {
  DepletionChain chain;
  ChainXmlDiagnostics diag;
  auto reactions = loadDepletionChainXmlString(chain, kChainXml, &diag);

  // Gd157 has no half_life -> registered but no decay data.
  EXPECT_GE(chain.indexOf({64, 157, 0}), 0);
  EXPECT_EQ(chain.decay({64, 157, 0}), nullptr);

  // The "(n,d)" reaction type is not modeled: dropped and counted.
  EXPECT_EQ(diag.unmodeledReactions, 1);
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

// Entries the reader cannot use are dropped and counted, never stored as a
// phantom nuclide.
TEST(ChainXml, UnusableEntriesAreDroppedAndCounted) {
  DepletionChain chain;
  ChainXmlDiagnostics diag;
  auto reactions = loadDepletionChainXmlString(chain, kDirtyChainXml, &diag);
  EXPECT_EQ(diag.unparsedNuclides, 1);
  EXPECT_EQ(diag.unparsedDecayTargets, 1);
  EXPECT_EQ(diag.unparsedYieldProducts, 1);
  EXPECT_EQ(diag.unmodeledReactions, 1);
  EXPECT_FALSE(diag.clean());
  EXPECT_TRUE(reactions.empty());
  EXPECT_EQ(chain.indexOf({0, 0, 0}), -1);

  const DecayData* i135 = chain.decay({53, 135, 0});
  ASSERT_NE(i135, nullptr);
  ASSERT_EQ(i135->modes.size(), 1u);  // the unparseable branch is gone
  EXPECT_EQ(i135->modes[0].daughter, (Zai{54, 135, 0}));
  EXPECT_DOUBLE_EQ(i135->modes[0].branching, 0.9);

  const FissionYields* y = chain.nearestYields({53, 135, 0}, 0.0253);
  ASSERT_NE(y, nullptr);
  EXPECT_EQ(y->products.size(), 1u);
}

TEST(ChainXml, LoadsFromAFile) {
  const TempFile f("cram_test_chain.xml", kChainXml);
  DepletionChain chain;
  auto reactions = loadDepletionChainXml(chain, f.path.string());
  EXPECT_GE(chain.indexOf({92, 235, 0}), 0);
  EXPECT_FALSE(reactions.empty());
}

TEST(ChainXml, MalformedAndMissingRootThrow) {
  DepletionChain chain;
  EXPECT_THROW(loadDepletionChainXml(chain, "/no/such/file_xyz.xml"), std::runtime_error);
  EXPECT_THROW(loadDepletionChainXmlString(chain, "<not_a_chain><foo/></not_a_chain>"),
               std::runtime_error);
  EXPECT_THROW(loadDepletionChainXmlString(chain, "<depletion_chain><unclosed>"),
               std::runtime_error);

  const TempFile bad("cram_test_bad_chain.xml", "<not_a_chain><foo/></not_a_chain>");
  EXPECT_THROW(loadDepletionChainXml(chain, bad.path.string()), std::runtime_error);
}

// A reaction whose product is not tracked (here U238 (n,p) -> "Nothing") must
// still remove the parent under irradiation.
TEST(ChainXml, UntrackedReactionProductStillRemovesParent) {
  DepletionChain chain;
  auto reactions = loadDepletionChainXmlString(chain, kChainXml);

  DepletionSystem sys(chain);
  // Give U238 only its (n,p)->Nothing channel a cross section.
  for (const auto& r : reactions) {
    if (r.parent == Zai{92, 238, 0} && r.type == ReactionType::NProton)
      sys.setReactions(r.parent, {reactionXs(r, 5.0)});
  }
  sys.setConstantFlux(1.0e16);

  Eigen::SparseMatrix<double> A = sys.assemble(Eigen::VectorXd::Zero(chain.size()));
  const int i = chain.indexOf({92, 238, 0});
  EXPECT_LT(A.coeff(i, i), 0.0);  // parent removed even though product untracked
}

// reactionXs() is how a chain's fission Q reaches the depletion system: a
// fission channel taken from the file satisfies constant-power normalization
// without the caller copying Q by hand, and a split channel's branching scales
// its cross section.
TEST(ChainXml, ReactionXsCarriesQAndBranching) {
  DepletionChain chain;
  auto reactions = loadDepletionChainXmlString(chain, kChainXml);

  DepletionSystem sys(chain);
  for (const auto& r : reactions) {
    if (r.parent == Zai{92, 235, 0} && r.type == ReactionType::Fission)
      sys.setReactions(r.parent, {reactionXs(r, 585.0)});
    if (r.parent == Zai{47, 109, 0} && r.target == Zai{47, 110, 1})
      sys.setReactions(r.parent, {reactionXs(r, 100.0)});
  }
  EXPECT_NO_THROW(sys.setConstantPower(1.0));

  bool sawMeta = false;
  for (const auto& [parent, rxns] : sys.reactions()) {
    if (parent != Zai{47, 109, 0})
      continue;
    ASSERT_EQ(rxns.size(), 1u);
    EXPECT_DOUBLE_EQ(rxns[0].sigma, 100.0 * 0.046);
    EXPECT_EQ(rxns[0].target, (Zai{47, 110, 1}));
    sawMeta = true;
  }
  EXPECT_TRUE(sawMeta);
}
