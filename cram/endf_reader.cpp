#include "cram/endf_reader.hpp"

#ifndef WITH_ENDFTK

namespace cram {
int loadDecayData(DepletionChain& /*chain*/, const std::string& /*path*/) {
  return 0;
}
int loadFissionYields(DepletionChain& /*chain*/, const std::string& /*path*/, bool /*cumulative*/) {
  return 0;
}
}  // namespace cram

#else
// =============================================================================
//  ENDFtk-backed implementation.
//
//  The accessor names below were verified against the ENDFtk source tree
//  (src/ENDFtk/section/8/{457,454,459}.hpp and their sub-records). Relevant
//  facts that the API pins down:
//    * section::Type<8,457>::halfLife() returns std::array<double,2>
//      = {value, uncertainty}  -> use [0].
//    * DecayMode::branchingRatio() returns a *range* {value, uncertainty}
//      -> take the first element.
//    * DecayMode::decayChain() is RTYP; finalIsomericState() is RFS.
//    * section::Type<8,457>::LISO() is the isomeric state of the PARENT.
//    * FissionYieldData::fissionProductIdentifiers()/isomericStates()/
//      fissionYieldValues() give the per-product ZAFP / FPS / Y columns.
//    * Columns are lazy ranges of double; materialise into a vector.
// =============================================================================
#include <ENDFtk.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

namespace cram {
namespace {

using namespace njoy::ENDFtk;

Zai zaiFromZA(int za, int liso) {
  return Zai{.z = za / 1000, .a = za % 1000, .i = liso};
}

// Materialise a lazy ENDFtk column range into a std::vector<double>.
template <typename Range>
std::vector<double> toVector(Range&& r) {
  std::vector<double> v;
  std::ranges::transform(std::forward<Range>(r), std::back_inserter(v),
                         [](auto x) { return static_cast<double>(x); });
  return v;
}

}  // namespace

int loadDecayData(DepletionChain& chain, const std::string& path) {
  auto tape = tree::fromFile(path);
  int count = 0;

  for (const auto& material : tape.materials()) {
    if (!material.hasSection(8, 457))
      continue;

    auto section = material.section(8, 457).parse<8, 457>();

    const int za = section.ZA();
    const int parentState = static_cast<int>(section.LISO());
    const Zai parent = zaiFromZA(za, parentState);

    DecayData d;
    d.halfLife = section.halfLife()[0];  // {value, uncertainty}

    if (!section.isStable()) {
      for (const auto& mode : section.decayModes().decayModes()) {
        DecayMode m;
        m.rtyp = mode.decayChain();                                  // RTYP
        m.finalState = static_cast<int>(mode.finalIsomericState());  // RFS

        // branchingRatio() is a {value, uncertainty} range.
        auto br = mode.branchingRatio();
        m.branching = static_cast<double>(*br.begin());

        bool fission = false;
        Zai daughter = applyDecay(parent, m.rtyp, m.finalState, fission);
        m.isFission = fission;
        d.modes.push_back(m);

        if (!fission)
          chain.add(daughter);  // register so production isn't dropped
      }
    }

    chain.setDecay(parent, std::move(d));
    ++count;
  }
  return count;
}

int loadFissionYields(DepletionChain& chain, const std::string& path, bool useCumulative) {
  const int mt = useCumulative ? 459 : 454;
  auto tape = tree::fromFile(path);
  int count = 0;

  for (const auto& material : tape.materials()) {
    if (!material.hasSection(8, mt))
      continue;

    // section::Type<8,454> and <8,459> share the same fission-yield layout.
    auto parseAndLoad = [&](const auto& section) {
      const Zai par = zaiFromZA(section.ZA(), 0);
      for (const auto& block : section.yields()) {
        FissionYields fy;
        fy.energy = block.incidentEnergy();

        auto zafp = toVector(block.fissionProductIdentifiers());
        auto fps = toVector(block.isomericStates());
        auto yld = toVector(block.fissionYieldValues());

        for (std::size_t k = 0; k < yld.size(); ++k) {
          Zai prod = zaiFromZA(static_cast<int>(std::lround(zafp[k])),
                               static_cast<int>(std::lround(fps[k])));
          fy.products.emplace_back(prod, yld[k]);
        }
        chain.addFissionYields(par, std::move(fy));
      }
    };

    if (useCumulative)
      parseAndLoad(material.section(8, 459).parse<8, 459>());
    else
      parseAndLoad(material.section(8, 454).parse<8, 454>());
    ++count;
  }
  return count;
}

}  // namespace cram
#endif  // WITH_ENDFTK
