#pragma once
//
// A seeded DepletionChain populated the way a real one is: nuclides carrying
// multi-mode decay data, plus a few fissioning parents with large independent
// yield tables.
//
// Matrix assembly is not a one-off cost. A model with many depletion regions
// rebuilds a burnup matrix per region, and does it again every burnup step, so
// decayMatrix() and the fission-source path run as often as the solver does.
// The benchmark suite previously measured only the solver, which left assembly
// unmeasured despite being on the same hot path.
//
#include <cmath>
#include <random>
#include <vector>

#include "cram/chain.hpp"

namespace cram_test {

using cram::DecayData;
using cram::DecayMode;
using cram::DepletionChain;
using cram::FissionYields;
using cram::Zai;

struct SyntheticDepletionChain {
  DepletionChain chain;
  std::vector<Zai> fissionParents;
};

// `n` nuclides with decay data, `fissionParents` of which also carry an
// independent fission-yield table of `productsPerParent` products -- the shape
// of an ENDF NFY table, which runs to ~1000 products for a major actinide.
inline SyntheticDepletionChain buildSyntheticDepletionChain(int n, int fissionParents = 8,
                                                            int productsPerParent = 800,
                                                            unsigned seed = 20260804u) {
  // Seeded as the other two fixtures are; see cyclic_chain.hpp for why the
  // single-value mt19937 constructor is avoided.
  std::seed_seq sequence{seed,        0x9e3779b9u, 0x85ebca6bu, 0xc2b2ae35u,
                         0x27d4eb2fu, 0x165667b1u, 0xd3a2646cu, 0xfd7046c5u};
  std::mt19937 rng(sequence);
  std::uniform_real_distribution<double> logHalfLife(1.0, 12.0);  // 10 s .. 1e12 s
  std::uniform_int_distribution<int> nModes(1, 3);
  std::uniform_int_distribution<int> rtypPick(0, 4);
  std::uniform_real_distribution<double> yield(1.0e-6, 6.0e-2);

  // Codes that move the daughter somewhere real: beta-, EC, alpha, IT, neutron.
  const double kRtyp[5] = {1.0, 2.0, 4.0, 3.0, 5.0};

  SyntheticDepletionChain out;

  // Lay the nuclides out over a plausible (Z, A) grid rather than a line, so
  // daughters land on genuinely different nuclides.
  std::vector<Zai> zais;
  zais.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const int z = 20 + (i % 76);   // Z = 20..95
    const int a = 40 + (i % 160);  // A = 40..199
    zais.push_back(Zai{.z = z, .a = a, .i = 0});
  }

  for (const Zai& parent : zais) {
    DecayData d;
    d.halfLife = std::pow(10.0, logHalfLife(rng));
    const int modes = nModes(rng);
    double remaining = 1.0;
    for (int k = 0; k < modes; ++k) {
      DecayMode m;
      m.rtyp = kRtyp[rtypPick(rng)];
      m.branching = (k == modes - 1) ? remaining : remaining / 2.0;
      remaining -= m.branching;
      m.finalState = 0;
      m.isFission = false;
      d.modes.push_back(m);
    }
    out.chain.setDecay(parent, std::move(d));
  }

  for (int f = 0; f < fissionParents && f < n; ++f) {
    const Zai parent = zais[static_cast<std::size_t>(f)];
    FissionYields fy;
    fy.energy = 0.0253;
    fy.products.reserve(static_cast<std::size_t>(productsPerParent));
    for (int p = 0; p < productsPerParent; ++p) {
      const Zai prod = zais[static_cast<std::size_t>((f * 7 + p * 3) % n)];
      fy.products.emplace_back(prod, yield(rng));
    }
    out.chain.addFissionYields(parent, std::move(fy));
    out.fissionParents.push_back(parent);
  }

  out.chain.close();  // register decay daughters so nothing is dropped
  return out;
}

}  // namespace cram_test
