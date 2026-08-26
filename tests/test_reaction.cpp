#include <gtest/gtest.h>

#include "cram/reaction.hpp"

using namespace cram;

TEST(Reaction, ProductTopology) {
  const Zai u238{92, 238, 0};
  EXPECT_EQ(reactionProduct(u238, ReactionType::NGamma), (Zai{92, 239, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::N2n), (Zai{92, 237, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::N3n), (Zai{92, 236, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::N4n), (Zai{92, 235, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::NAlpha), (Zai{90, 235, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::NProton), (Zai{91, 238, 0}));
  EXPECT_EQ(reactionProduct(u238, ReactionType::Fission), (Zai{92, 238, 0}));  // unchanged
}

// The product is always the ground state, whatever the parent's level.
TEST(Reaction, ProductIsGroundState) {
  const Zai am242m{95, 242, 1};
  EXPECT_EQ(reactionProduct(am242m, ReactionType::NGamma), (Zai{95, 243, 0}));
}

TEST(Reaction, NamesRoundTrip) {
  for (ReactionType t :
       {ReactionType::Fission, ReactionType::NGamma, ReactionType::N2n, ReactionType::N3n,
        ReactionType::N4n, ReactionType::NAlpha, ReactionType::NProton}) {
    const auto back = reactionTypeFromName(reactionName(t));
    ASSERT_TRUE(back.has_value()) << reactionName(t);
    EXPECT_EQ(*back, t);
  }
}

TEST(Reaction, UnknownNameIsEmpty) {
  EXPECT_FALSE(reactionTypeFromName("(n,d)").has_value());
  EXPECT_FALSE(reactionTypeFromName("").has_value());
  EXPECT_FALSE(reactionTypeFromName("Fission").has_value());  // case-sensitive
}
