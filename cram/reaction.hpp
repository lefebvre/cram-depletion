#pragma once
//
// Neutron-induced transmutation reactions: the channel types a depletion
// calculation tracks, their conventional ground-state products, and the
// reaction names OpenMC's depletion chain uses for them.
//
#include <optional>
#include <string_view>

#include "cram/nuclide.hpp"

namespace cram {

enum class ReactionType { Fission, NGamma, N2n, N3n, N4n, NAlpha, NProton };

// The ground-state product of `type` on `parent`: (n,gamma) A+1, (n,2n) A-1,
// (n,3n) A-2, (n,4n) A-3, (n,alpha) Z-2/A-3, (n,p) Z-1. Fission has no single
// product and returns `parent` unchanged; callers branch on the type first.
// Metastable products are not derived here -- a chain that tracks them names
// the target explicitly (see ReactionXS::target).
Zai reactionProduct(const Zai& parent, ReactionType type);

// OpenMC depletion-chain reaction names: "fission", "(n,gamma)", "(n,2n)",
// "(n,3n)", "(n,4n)", "(n,a)", "(n,p)". Empty for a name that is not one of
// the channels tracked here.
std::optional<ReactionType> reactionTypeFromName(std::string_view name);
const char* reactionName(ReactionType type);

// One neutron reaction channel of a chain's topology, without its cross
// section: which parent, which channel, which product (if any is tracked), the
// reaction Q value, and the fraction of the channel that goes to this product
// (a chain that tracks isomers splits (n,gamma) into a ground-state and a
// metastable entry whose branchings sum to 1). This is what the OpenMC chain
// reader returns and what a cross-section source is matched against.
//
// No default member initializers: a channel with a forgotten parent or type
// is a phantom entry, a forgotten Q makes constant-power normalization wrong
// for a fission channel, and a forgotten branching double-counts a split
// channel. `target` is std::nullopt when the product is not tracked.
struct ChainReaction {
  Zai parent;
  ReactionType type;
  std::optional<Zai> target;
  double q;          // reaction Q value [eV]
  double branching;  // fraction of the channel going to `target`
};

}  // namespace cram
