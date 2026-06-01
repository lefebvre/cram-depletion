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
#include <string>
#include <vector>

#include "cram/chain.hpp"
#include "cram/deplete.hpp"

namespace cram {

// One neutron reaction channel parsed from the chain (cross section excluded).
struct ChainReaction {
  Zai parent;
  ReactionType type = ReactionType::NGamma;
  Zai target{};  // {0,0,0} when the reaction has no tracked product
  bool hasTarget = false;
  double q = 0.0;          // reaction Q value [eV]
  double branching = 1.0;  // branch fraction (e.g. (n,gamma) to ground vs metastable)
};

// Parse a GND-style nuclide name ("U235", "Am242_m1", "Th231") into a Zai.
// Returns false if the element symbol or mass cannot be recognized.
bool parseNuclideName(const std::string& name, Zai& out);

// Element symbol for proton number z (1..118), or empty string if out of range.
std::string elementSymbol(int z);

// Parse a depletion_chain XML file into `chain` (registers nuclides, sets decay
// data and fission yields) and return the neutron reaction topology. Throws
// std::runtime_error if the file cannot be opened or parsed.
std::vector<ChainReaction> loadDepletionChainXml(DepletionChain& chain, const std::string& path);

}  // namespace cram
