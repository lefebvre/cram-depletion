#include "cram/reaction.hpp"

namespace cram {

Zai reactionProduct(const Zai& parent, ReactionType type) {
  Zai p{.z = parent.z, .a = parent.a, .i = 0};
  switch (type) {
    case ReactionType::NGamma:  // (n,gamma): A -> A+1
      p.a += 1;
      break;
    case ReactionType::N2n:  // (n,2n): A -> A-1
      p.a -= 1;
      break;
    case ReactionType::N3n:  // (n,3n): A -> A-2
      p.a -= 2;
      break;
    case ReactionType::N4n:  // (n,4n): A -> A-3
      p.a -= 3;
      break;
    case ReactionType::NAlpha:  // (n,alpha): Z-2, A-3
      p.z -= 2;
      p.a -= 3;
      break;
    case ReactionType::NProton:  // (n,p): Z-1, A unchanged
      p.z -= 1;
      break;
    case ReactionType::Fission:  // no single product
      break;
  }
  return p;
}

std::optional<ReactionType> reactionTypeFromName(std::string_view name) {
  if (name == "fission")
    return ReactionType::Fission;
  if (name == "(n,gamma)")
    return ReactionType::NGamma;
  if (name == "(n,2n)")
    return ReactionType::N2n;
  if (name == "(n,3n)")
    return ReactionType::N3n;
  if (name == "(n,4n)")
    return ReactionType::N4n;
  if (name == "(n,a)")
    return ReactionType::NAlpha;
  if (name == "(n,p)")
    return ReactionType::NProton;
  return std::nullopt;
}

const char* reactionName(ReactionType type) {
  switch (type) {
    case ReactionType::Fission:
      return "fission";
    case ReactionType::NGamma:
      return "(n,gamma)";
    case ReactionType::N2n:
      return "(n,2n)";
    case ReactionType::N3n:
      return "(n,3n)";
    case ReactionType::N4n:
      return "(n,4n)";
    case ReactionType::NAlpha:
      return "(n,a)";
    case ReactionType::NProton:
      return "(n,p)";
  }
  return "";
}

}  // namespace cram
