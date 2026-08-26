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
// the file but dropped: a nuclide whose name did not parse (its whole element
// is skipped), a decay target that did not parse (that mode is skipped -- its
// branching is lost, not rerouted), a yield product that did not parse, or a
// reaction of a type this library does not model. Handed back through the
// out-parameter of the loaders when the caller asks; written to stderr as a
// warning otherwise, so it can never pass unnoticed either way.
struct ChainXmlDiagnostics {
  int unparsedNuclides = 0;
  int unparsedDecayTargets = 0;
  int unparsedYieldProducts = 0;
  int unmodeledReactions = 0;

  bool clean() const {
    return unparsedNuclides == 0 && unparsedDecayTargets == 0 && unparsedYieldProducts == 0 &&
           unmodeledReactions == 0;
  }
};

// Parse a depletion_chain XML file into `chain` (registers nuclides, sets decay
// data and fission yields) and return the neutron reaction topology. Every
// decay mode carries its explicit target as DecayMode::daughter, so the chain
// honors the file's topology even where the RTYP rules would differ. Throws
// std::runtime_error if the file cannot be opened or parsed, or if the root
// element is not <depletion_chain>.
std::vector<ChainReaction> loadDepletionChainXml(DepletionChain& chain, const std::string& path,
                                                 ChainXmlDiagnostics* diagnostics = nullptr);

// Same, from an in-memory document.
std::vector<ChainReaction> loadDepletionChainXmlString(DepletionChain& chain, std::string_view xml,
                                                       ChainXmlDiagnostics* diagnostics = nullptr);

}  // namespace cram
