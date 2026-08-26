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

}  // namespace cram
