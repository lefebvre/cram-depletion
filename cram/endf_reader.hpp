#pragma once
//
// Ingestion of ENDF/B-VIII data into a DepletionChain using njoy/ENDFtk.
//
// Compiled only with CRAM_WITH_ENDFTK (it needs njoy/ENDFtk). Without it the
// loaders throw std::runtime_error, as the OpenMC chain reader's do, so a caller
// built without the reader finds out at the call rather than by a chain that
// silently stayed empty. The ENDFtk accessor names used in the implementation
// have been verified against the ENDFtk source tree and a real build (see
// tests/integration).
//
// ENDF layout used here:
//   * Decay sublibrary:        MF8 / MT457   (half-life, decay modes, BRs)
//   * Neutron fission yields:  MF8 / MT454   (INDEPENDENT yields)  <-- use these
//                              MF8 / MT459   (cumulative yields)   <-- do NOT mix
//                                            with an explicit decay chain.
//
#include <string>

#include "cram/chain.hpp"

namespace cram {

// Read every MF8/MT457 section from an ENDF decay-data tape into the chain.
// Returns the number of nuclides whose decay data was loaded.
int loadDecayData(DepletionChain& chain, const std::string& endfPath);

// Read MF8/MT454 (independent) fission yields from an ENDF NFY tape.
// Pass useCumulative=true to read MT459 instead (only if you are NOT modeling
// the decay chain explicitly).
int loadFissionYields(DepletionChain& chain, const std::string& endfPath,
                      bool useCumulative = false);

}  // namespace cram
