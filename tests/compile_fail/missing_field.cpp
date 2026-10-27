// Built only by the CompileFail tests in tests/CMakeLists.txt, never by a
// normal build.
//
// Several structs leave their required members without a default member
// initializer, so that -Wmissing-field-initializers turns a forgotten field
// into a diagnostic instead of a plausible-looking zero (see DecayMode in
// chain.hpp). That guarantee lives in the compiler, not the type system, so
// it is tested: with CRAM_OMIT unset every initialization below names each
// required field and must compile cleanly; with CRAM_OMIT=<n>, initialization
// <n> leaves one required field out and the build must fail on exactly that
// warning.
#include <optional>

#include "cram/chain.hpp"
#include "cram/deplete.hpp"
#include "cram/nuclide.hpp"
#include "cram/reaction.hpp"

#ifndef CRAM_OMIT
#define CRAM_OMIT 0
#endif

int main() {
  using namespace cram;
#if CRAM_OMIT == 1
  const Zai z{.z = 92};  // no mass number
#else
  const Zai z{.z = 92, .a = 235};
#endif
#if CRAM_OMIT == 2
  const DecayMode m{.rtyp = 1.0, .branching = 1.0, .finalState = 0};  // no isFission
#else
  const DecayMode m{.rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false};
#endif
#if CRAM_OMIT == 3
  const DecayData d{.gammaEnergyPerDecay = 0.0};  // no halfLife
#else
  const DecayData d{.halfLife = 1.0};
#endif
#if CRAM_OMIT == 4
  const FissionYields y{.products = {}};  // no energy
#else
  const FissionYields y{.energy = 0.0253};
#endif
#if CRAM_OMIT == 5
  const ReactionXS r{.type = ReactionType::NGamma, .target = std::nullopt};  // no sigma
#else
  const ReactionXS r{.type = ReactionType::NGamma, .target = std::nullopt, .sigma = 1.0};
#endif
#if CRAM_OMIT == 6
  const ChainReaction c{.parent = z,
                        .type = ReactionType::Fission,
                        .target = std::nullopt,
                        .q = 2.0e8};  // no branching
#else
  const ChainReaction c{.parent = z,
                        .type = ReactionType::Fission,
                        .target = std::nullopt,
                        .q = 2.0e8,
                        .branching = 1.0};
#endif
  return (m.isFission || d.halfLife < 0.0 || y.energy < 0.0 || r.sigma < 0.0 || c.branching < 0.0 ||
          !z.representable())
             ? 1
             : 0;
}
