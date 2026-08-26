#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>

#include "cram/chain.hpp"
#include "cram/chain_xml.hpp"
#include "cram/cram.hpp"
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
  // explicit daughter, and an SF branch back to U235 itself -- the target
  // OpenMC writes on an sf mode, honored like any other rather than routed to
  // the yield tables.
  const DecayData* u235 = chain.decay({92, 235, 0});
  ASSERT_NE(u235, nullptr);
  EXPECT_NEAR(u235->halfLife, 2.22102e16, 1e9);
  const auto alpha = modeTo(*u235, Zai{90, 231, 0});
  ASSERT_TRUE(alpha.has_value());
  EXPECT_DOUBLE_EQ(alpha->branching, 0.99999999993);
  EXPECT_FALSE(alpha->isFission);
  const auto sf = modeTo(*u235, Zai{92, 235, 0});
  ASSERT_TRUE(sf.has_value());
  EXPECT_DOUBLE_EQ(sf->branching, 7.0e-11);
  EXPECT_FALSE(sf->isFission);
  EXPECT_EQ(u235->modes.size(), 2u);

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

// OpenMC delegates the yields of a nuclide with no measured ones of its own to
// another nuclide, <neutron_fission_yields parent="U235"/>, and does so for a
// good fraction of the fissionable nuclides in a real chain. Read past, those
// nuclides fission with no products and are never burned at all.
TEST(ChainXml, DelegatedFissionYieldsAreResolved) {
  const char* xml = R"XML(<depletion_chain>
  <nuclide name="U235" half_life="2.22102e16" decay_modes="0" reactions="1">
    <reaction type="fission" Q="193405400.0"/>
    <neutron_fission_yields>
      <energies>0.0253</energies>
      <fission_yields energy="0.0253">
        <products>Xe135 I135</products>
        <data>0.00254 0.0292</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="Np236" half_life="4.8e12" reactions="1">
    <reaction type="fission" Q="190000000.0"/>
    <neutron_fission_yields parent="U235"/>
  </nuclide>
  <nuclide name="Pu237" reactions="1">
    <reaction type="fission" Q="190000000.0"/>
    <neutron_fission_yields parent="Np236"/>
  </nuclide>
  <nuclide name="Th230" reactions="1">
    <reaction type="fission" Q="180000000.0"/>
    <neutron_fission_yields parent="Cf252"/>
  </nuclide>
  <nuclide name="Xe135"/>
  <nuclide name="I135"/>
</depletion_chain>)XML";

  DepletionChain chain;
  ChainXmlDiagnostics diag;
  const auto reactions = loadDepletionChainXmlString(chain, xml, &diag);

  const FissionYields* source = chain.nearestYields({92, 235, 0}, 0.0253);
  ASSERT_NE(source, nullptr);
  // Np236 delegates to U235; Pu237 delegates to Np236, which delegates on.
  for (const Zai& z : {Zai{93, 236, 0}, Zai{94, 237, 0}}) {
    const FissionYields* y = chain.nearestYields(z, 0.0253);
    ASSERT_NE(y, nullptr) << z.str();
    EXPECT_DOUBLE_EQ(y->energy, source->energy) << z.str();
    ASSERT_EQ(y->products.size(), source->products.size()) << z.str();
    for (std::size_t k = 0; k < y->products.size(); ++k) {
      EXPECT_EQ(y->products[k].first, source->products[k].first) << z.str();
      EXPECT_DOUBLE_EQ(y->products[k].second, source->products[k].second) << z.str();
    }
  }

  // A delegation to a nuclide the file does not carry cannot be resolved, and
  // is counted rather than leaving a silently yield-less nuclide behind.
  EXPECT_EQ(chain.nearestYields({90, 230, 0}, 0.0253), nullptr);
  EXPECT_EQ(diag.unresolvedYieldDelegations, 1);
  EXPECT_FALSE(diag.clean());

  // The delegated table is what the burnup matrix uses: Np236 fission produces
  // the fission products, not nothing.
  DepletionSystem sys(chain);
  for (const auto& r : reactions) {
    if (r.parent == Zai{93, 236, 0} && r.type == ReactionType::Fission)
      sys.setReactions(r.parent, {reactionXs(r, 10.0)});
  }
  sys.setConstantFlux(1.0e14);
  const Eigen::SparseMatrix<double> A = sys.assemble(Eigen::VectorXd::Zero(chain.size()));
  const int np = chain.indexOf({93, 236, 0});
  EXPECT_LT(A.coeff(np, np), 0.0);
  EXPECT_GT(A.coeff(chain.indexOf({54, 135, 0}), np), 0.0);
  EXPECT_GT(A.coeff(chain.indexOf({53, 135, 0}), np), 0.0);
}

// An OpenMC chain names the parent itself as the target of a spontaneous-
// fission mode, so the branch removes and restores the same atoms. Reading it
// as a fission mode instead removes the parent at branching*lambda and
// substitutes the nearest NEUTRON-induced yield set for its products -- and
// where the parent has no yield table at all, the atoms simply vanish, with
// decayMatrix()'s dropped-daughter counter none the wiser.
TEST(ChainXml, SpontaneousFissionSelfTargetLeavesTheParentAlone) {
  const char* xml = R"XML(<depletion_chain>
  <nuclide name="U232" half_life="2.17e9" decay_modes="1" reactions="0">
    <decay type="sf" target="U232" branching_ratio="1.0"/>
  </nuclide>
</depletion_chain>)XML";

  DepletionChain chain;
  loadDepletionChainXmlString(chain, xml);
  EXPECT_EQ(chain.close(), 0);
  int dropped = -1;
  const Eigen::SparseMatrix<double> M = chain.decayMatrix(&dropped);
  EXPECT_EQ(dropped, 0);
  const int i = chain.indexOf({92, 232, 0});
  ASSERT_GE(i, 0);
  EXPECT_DOUBLE_EQ(M.coeff(i, i), 0.0) << "removal and self-production cancel";

  Eigen::VectorXd n0(chain.size());
  n0.setZero();
  n0(i) = 1.0;
  EXPECT_NEAR(cramSolve(M, n0, 1.0e10)(i), 1.0, 1e-12) << "and nothing is lost over any span";
}

// OpenMC omits the target of a decay mode whenever the daughter falls outside
// the chain, and spells an untracked reaction product "Nothing". Neither is an
// entry the reader failed on, so neither may make clean() false: a caller that
// passes no diagnostics pointer would get a warning about a perfectly good
// file, and one that checks the counters would see a failure that is not there.
TEST(ChainXml, AbsentDecayTargetIsNotAParseFailure) {
  const char* xml = R"XML(<depletion_chain>
  <nuclide name="Cm244" half_life="5.715e8" decay_modes="2" reactions="1">
    <decay type="alpha" target="Pu240" branching_ratio="0.9999"/>
    <decay type="beta-" branching_ratio="0.0001"/>
    <reaction type="(n,gamma)" Q="6800000.0" target="Nothing"/>
  </nuclide>
  <nuclide name="Pu240"/>
</depletion_chain>)XML";

  DepletionChain chain;
  ChainXmlDiagnostics diag;
  loadDepletionChainXmlString(chain, xml, &diag);
  EXPECT_EQ(diag.unparsedDecayTargets, 0);
  EXPECT_EQ(diag.unparsedReactionTargets, 0);
  EXPECT_TRUE(diag.clean());

  // The targetless branch is dropped, but the parent still leaves at the full
  // decay constant -- only its 0.01% of production is missing.
  const DecayData* cm = chain.decay({96, 244, 0});
  ASSERT_NE(cm, nullptr);
  ASSERT_EQ(cm->modes.size(), 1u);
  const Eigen::SparseMatrix<double> M = chain.decayMatrix();
  const int i = chain.indexOf({96, 244, 0});
  EXPECT_DOUBLE_EQ(M.coeff(i, i), -cm->decayConstant);
  EXPECT_NEAR(M.coeff(chain.indexOf({94, 240, 0}), i), 0.9999 * cm->decayConstant,
              1e-12 * cm->decayConstant);
}

// A reaction target the file names but never declares as its own <nuclide> is
// registered anyway. Left out of the chain it is indistinguishable from the
// deliberate "Nothing": assemble() consumes the parent and produces nothing,
// and no counter, warning or close() notices -- the one hole in the reader's
// promise that a dropped entry can never pass unnoticed. A target that cannot
// be read at all is what the counter is for.
TEST(ChainXml, ReactionTargetsAreRegisteredAndUnreadableOnesCounted) {
  const char* xml = R"XML(<depletion_chain>
  <nuclide name="Ag109" reactions="3">
    <reaction type="(n,gamma)" Q="6809000.0" target="Ag110" branching_ratio="0.954"/>
    <reaction type="(n,gamma)" Q="6809000.0" target="Ag110_m1" branching_ratio="0.046"/>
    <reaction type="(n,2n)" Q="-9000000.0" target="Xx108"/>
  </nuclide>
</depletion_chain>)XML";

  DepletionChain chain;
  ChainXmlDiagnostics diag;
  const auto reactions = loadDepletionChainXmlString(chain, xml, &diag);
  EXPECT_GE(chain.indexOf({47, 110, 0}), 0);
  EXPECT_GE(chain.indexOf({47, 110, 1}), 0);
  EXPECT_EQ(chain.decay({47, 110, 0}), nullptr) << "registered bare, carrying no data";
  EXPECT_EQ(diag.unparsedReactionTargets, 1);  // Xx108
  EXPECT_FALSE(diag.clean());

  DepletionSystem sys(chain);
  for (const auto& r : reactions) {
    if (r.type == ReactionType::NGamma)
      sys.setReactions(r.parent, {reactionXs(r, 91.0)});
  }
  sys.setConstantFlux(1.0e14);
  const Eigen::SparseMatrix<double> A = sys.assemble(Eigen::VectorXd::Zero(chain.size()));
  const int ag = chain.indexOf({47, 109, 0});
  // Whatever leaves Ag109 arrives at the metastable branch: nothing is lost.
  double colSum = 0.0;
  for (int r = 0; r < A.rows(); ++r)
    colSum += A.coeff(r, ag);
  EXPECT_NEAR(colSum, 0.0, 1e-20);
  EXPECT_GT(A.coeff(chain.indexOf({47, 110, 1}), ag), 0.0);
}
