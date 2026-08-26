#pragma once
//
// Reader for OpenMC's depletion_chain XML (the file openmc.deplete.Chain
// exports / consumes). It carries the decay data, fission yields, and the
// neutron reaction topology of a depletion chain, but NOT cross sections --
// those are supplied separately (from transport, or a micro-xs file). This is
// the bridge that lets this engine consume a real CASL/VERA chain.
//
// The schema (see Yu & Forget 2022, Fig. 1):
//   <depletion_chain>
//     <nuclide name="U235" half_life="..." decay_modes="2" reactions="5">
//       <decay type="alpha" target="Th231" branching_ratio="0.99..."/>
//       <decay type="sf" target="U235" branching_ratio="7.0e-11"/>
//       <reaction type="(n,gamma)" Q="6545200.0" target="U236"/>
//       <reaction type="fission" Q="193405400.0"/>
//       <neutron_fission_yields>
//         <energies>0.0253</energies>
//         <fission_yields energy="0.0253">
//           <products>Xe135 I135 ...</products>
//           <data>0.00254 0.0292 ...</data>
//         </fission_yields>
//       </neutron_fission_yields>
//     </nuclide>
//   </depletion_chain>
//
// Two spellings in that schema are easy to read past, and both change the
// matrix:
//
//   * A spontaneous-fission <decay type="sf"> carries a target like any other
//     mode, and OpenMC writes the parent itself there -- the mode removes and
//     restores the same atoms, so an `sf` branch of a chain whose fission
//     products are not tracked as such changes nothing. The target is honored
//     as the mode's daughter. Only a targetless `sf` mode falls back to the
//     SFY table (yields at energy 0), which is how a chain assembled in this
//     library rather than read from OpenMC expresses spontaneous fission.
//
//   * <neutron_fission_yields parent="U235"/> delegates a nuclide's yields to
//     another nuclide's tables instead of repeating them; OpenMC writes it for
//     every fissionable nuclide with no measured yields of its own. The
//     delegated tables are resolved (through further delegations) and stored
//     for the delegating nuclide, so its fission channel has products.
//
// Compiled only with CRAM_WITH_CHAIN_XML (it needs pugixml). Without it the
// loaders throw std::runtime_error, so a caller built without the reader finds
// out at the call rather than by a chain that silently stayed empty.
//
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "cram/chain.hpp"
#include "cram/reaction.hpp"

namespace cram {

// Parse a GND-style nuclide name ("U235", "Am242_m1", "Th231") into a Zai.
// Empty if the element symbol or mass cannot be recognized.
std::optional<Zai> parseNuclideName(std::string_view name);

// Element symbol for proton number z (1..118), or empty string if out of range.
std::string elementSymbol(int z);

// What the reader could not use. Each count is an entry that was present in
// the file but dropped:
//
//   * a nuclide whose name did not parse -- its whole element is skipped;
//   * a decay target that did not parse -- that mode is skipped, its branching
//     lost, not rerouted (an ABSENT target, or "Nothing", is not a failure:
//     OpenMC writes it whenever the daughter falls outside the chain, and the
//     parent is still removed at the full decay constant);
//   * a fission-yield product that did not parse;
//   * a reaction of a type this library does not model;
//   * a <neutron_fission_yields parent="..."/> whose target has no tables to
//     delegate -- that nuclide is left with no yields, so its fission channel
//     consumes it and produces nothing;
//   * a reaction target that did not parse -- the channel is kept and still
//     consumes its parent, but its product is lost (again, an absent target or
//     "Nothing" is the file's deliberate statement that the product is not
//     tracked, not a failure).
//
// Handed back through the out-parameter of the loaders when the caller asks;
// written to stderr as a warning otherwise, so it can never pass unnoticed
// either way.
struct ChainXmlDiagnostics {
  int unparsedNuclides = 0;
  int unparsedDecayTargets = 0;
  int unparsedYieldProducts = 0;
  int unmodeledReactions = 0;
  int unresolvedYieldDelegations = 0;
  int unparsedReactionTargets = 0;

  bool clean() const {
    return unparsedNuclides == 0 && unparsedDecayTargets == 0 && unparsedYieldProducts == 0 &&
           unmodeledReactions == 0 && unresolvedYieldDelegations == 0 &&
           unparsedReactionTargets == 0;
  }
};

// Parse a depletion_chain XML file into `chain` (registers nuclides, sets decay
// data and fission yields) and return the neutron reaction topology. Every
// decay mode carries its explicit target as DecayMode::daughter, so the chain
// honors the file's topology even where the RTYP rules would differ.
//
// A reaction target the file names but never declares as its own <nuclide> is
// registered too, as a bare terminator carrying no data -- the same thing
// DepletionChain::close() does for an unregistered decay daughter, and for the
// same reason: an unregistered target is one matrix assembly consumes the
// parent for and produces nothing from, indistinguishable from the deliberate
// "Nothing". A target that is absent or spelled "Nothing" stays untracked,
// which is what the file is asking for.
//
// Throws std::runtime_error if the file cannot be opened or parsed, or if the
// root element is not <depletion_chain>.
std::vector<ChainReaction> loadDepletionChainXml(DepletionChain& chain, const std::string& path,
                                                 ChainXmlDiagnostics* diagnostics = nullptr);

// Same, from an in-memory document.
std::vector<ChainReaction> loadDepletionChainXmlString(DepletionChain& chain, std::string_view xml,
                                                       ChainXmlDiagnostics* diagnostics = nullptr);

}  // namespace cram
