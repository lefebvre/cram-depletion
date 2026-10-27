// Integration test: exercises the real ENDFtk-backed reader against a known
// MT457 decay section (Am-242m) and a real MT454 fission-yield file (Pu-239).
// Built only with CRAM_WITH_ENDFTK.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "cram/chain.hpp"
#include "cram/endf_reader.hpp"

using namespace cram;

// Am-242m1: ZA=95242, LISO=1, half-life 4.449622e9 s.
// Three decay branches: alpha (RTYP 4), IT (RTYP 3), spontaneous fission
// (RTYP 6); IT dominates (BR ~0.9954).
TEST(EndfReaderIntegration, Am242mDecay) {
  DepletionChain chain;
  int n = loadDecayData(chain, "data/decay_am242m.endf");
  ASSERT_EQ(n, 1) << "expected one MT457 section";

  const Zai am242m{95, 242, 1};
  const DecayData* d = chain.decay(am242m);
  ASSERT_NE(d, nullptr) << "Am-242m decay data should be present";
  EXPECT_NEAR(d->halfLife, 4.449622e9, 1.0);
  ASSERT_EQ(d->modes.size(), 3u) << "Am-242m has three decay branches";

  const double lam = kLn2 / 4.449622e9;
  EXPECT_NEAR(d->decayConstant, lam, lam * 1e-12);

  bool sawFission = false;
  bool sawAlpha = false;
  double brSum = 0.0;
  for (const auto& m : d->modes) {
    brSum += m.branching;
    if (m.isFission)
      sawFission = true;
    if (std::abs(m.rtyp - 4.0) < 1e-9)
      sawAlpha = true;
  }
  EXPECT_TRUE(sawFission) << "spontaneous-fission branch should be flagged";
  EXPECT_TRUE(sawAlpha) << "alpha branch should be present";
  EXPECT_GT(brSum, 0.9);
  EXPECT_LT(brSum, 1.1);
}

TEST(EndfReaderIntegration, Pu239ThermalFissionYields) {
  DepletionChain chain;
  int n = loadFissionYields(chain, "data/nfy-Pu239.endf", /*useCumulative=*/false);
  ASSERT_EQ(n, 1) << "expected one MT454 material";

  const Zai pu239{94, 239, 0};
  const FissionYields* thermal = chain.nearestYields(pu239, 0.0253);
  ASSERT_NE(thermal, nullptr) << "thermal yields present";
  EXPECT_FALSE(thermal->products.empty());

  // Sum of independent yields per fission should be ~2 (two fragments).
  double ySum = 0.0;
  for (const auto& [zai, y] : thermal->products)
    ySum += y;
  EXPECT_GT(ySum, 1.8);
  EXPECT_LT(ySum, 2.2);

  // The matrix builder must accept a fission source built from these yields.
  chain.add(pu239);
  std::vector<Eigen::Triplet<double>> t;
  chain.addFissionSource(t, pu239, /*fissionRate=*/1e-3, /*energy=*/0.0253);
  auto A = chain.finalize(t);
  EXPECT_GT(A.nonZeros(), 1);
}

// The parent's isomeric state lives in MF1/MT451 (LISO), not in MT454's HEAD
// record. Rewrite LISO in a copy of the Pu-239 tape to 1, as an Am-242m1 tape
// carries it, and the yields must land on the isomer, not on the ground state.
TEST(EndfReaderIntegration, IsomerParentYieldsAreFiledUnderTheIsomer) {
  std::ifstream in("data/nfy-Pu239.endf");
  ASSERT_TRUE(in.good());
  std::vector<std::string> lines;
  for (std::string line; std::getline(in, line);)
    lines.push_back(line);
  // Line 3 is MT451's second record: ELIS, STA, LIS, LISO, 0, NFOR in 11-column
  // fields, so LISO occupies columns 34-44.
  ASSERT_GT(lines.size(), 3u);
  ASSERT_GE(lines[2].size(), 75u);
  ASSERT_EQ(lines[2].substr(66, 9), "9437 1451") << "MT451's second record moved";
  lines[2].replace(33, 11, "          1");

  const std::filesystem::path copy =
      std::filesystem::temp_directory_path() / "cram_nfy_Pu239m1.endf";
  {
    std::ofstream out(copy);
    for (const std::string& line : lines)
      out << line << '\n';
  }
  DepletionChain chain;
  const int n = loadFissionYields(chain, copy.string());
  std::filesystem::remove(copy);
  ASSERT_EQ(n, 1);
  EXPECT_NE(chain.nearestYields({94, 239, 1}, 0.0253), nullptr);
  EXPECT_EQ(chain.nearestYields({94, 239, 0}, 0.0253), nullptr);
}
