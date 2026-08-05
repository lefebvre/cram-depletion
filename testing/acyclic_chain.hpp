#pragma once
//
// The decay-only counterpart to cyclic_chain.hpp: a seeded synthetic burnup
// matrix whose transitions all run forward, so the nuclides are already in a
// topological order and the matrix is acyclic.
//
// This is not a degenerate case -- it is what a pure-decay problem looks like.
// Radioactive decay always runs downhill, so a decay chain can always be
// topologically ordered; cycles enter a burnup matrix only once neutron
// reactions do, because (n,gamma) walks the chain up while decay walks it back
// down. DepletionChain::decayMatrix() produces exactly this shape, and it is
// what cram-apps/deplete.cpp and the whole ENDF MT457 path operate on.
//
// SparseLU exploits it automatically: an acyclic matrix has an elimination
// order with little or no fill-in and COLAMD finds it, so decay-only problems
// are roughly twice as fast as burnup problems of the same size with no solver
// option involved. Both regimes are benchmarked so that a change which helps
// one is not assumed to help the other.
//
// Deliberately kept separate from cyclic_chain.hpp rather than folded into it
// with a flag: that generator's exact output is frozen in
// tests/cram_golden_data.hpp, so it must not be disturbed. Everything here
// matches it -- same decay-constant range, same branching normalization, same
// initial inventory -- except the choice of daughter index, so the two differ
// only in topology and stay directly comparable.
//
#include <Eigen/SparseCore>
#include <cmath>
#include <random>
#include <vector>

namespace cram_test {

struct AcyclicChain {
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd n0;
};

inline AcyclicChain buildAcyclicChain(int n, unsigned seed = 20260804u) {
  // Seeded exactly as cyclic_chain.hpp is, and for the same reason: mt19937's
  // single-value constructor reaches only 2^32 of its 2^19937 initial states
  // and leaves them structured rather than mixed. See that file for the full
  // reasoning. Matching the two matters here beyond tidiness -- these fixtures
  // are compared against each other in the benchmarks, so they should differ in
  // topology and nothing else.
  std::seed_seq sequence{seed,        0x9e3779b9u, 0x85ebca6bu, 0xc2b2ae35u,
                         0x27d4eb2fu, 0x165667b1u, 0xd3a2646cu, 0xfd7046c5u};
  std::mt19937 rng(sequence);
  std::uniform_real_distribution<double> logLambda(-9.0, -1.0);  // 1e-9 .. 1e-1 /s
  std::uniform_int_distribution<int> nModes(1, 3);
  std::uniform_int_distribution<int> hop(1, 4);
  std::uniform_real_distribution<double> branch(0.15, 1.0);
  std::uniform_int_distribution<int> stableRoll(0, 19);  // ~5% stable

  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n) * 5);

  for (int i = 0; i < n; ++i) {
    if (i != 0 && stableRoll(rng) == 0)
      continue;  // stable sink: no diagonal, no production

    const double lambda = std::pow(10.0, logLambda(rng));
    const int modes = nModes(rng);

    // Every daughter sits at a higher index than its parent, so no cycle can
    // form. The last nuclide is a terminator and absorbs anything past the end.
    std::vector<int> targets;
    targets.reserve(static_cast<std::size_t>(modes));
    for (int k = 0; k < modes; ++k) {
      int j = i + hop(rng);
      if (j >= n)
        j = n - 1;
      if (j == i)
        j = (i + 1 < n) ? i + 1 : i;
      targets.push_back(j);
    }

    std::vector<double> br(static_cast<std::size_t>(modes));
    double drawn = 0.0;
    for (int k = 0; k < modes - 1; ++k) {
      br[static_cast<std::size_t>(k)] = branch(rng) / modes;
      drawn += br[static_cast<std::size_t>(k)];
    }
    if (drawn >= 1.0) {
      for (int k = 0; k < modes; ++k)
        br[static_cast<std::size_t>(k)] = 1.0 / modes;
    } else {
      br[static_cast<std::size_t>(modes - 1)] = 1.0 - drawn;
    }

    t.emplace_back(i, i, -lambda);
    for (int k = 0; k < modes; ++k)
      t.emplace_back(targets[static_cast<std::size_t>(k)], i,
                     lambda * br[static_cast<std::size_t>(k)]);
  }

  AcyclicChain c;
  c.A.resize(n, n);
  c.A.setFromTriplets(t.begin(), t.end());  // sums duplicates
  c.A.makeCompressed();

  std::uniform_real_distribution<double> logN0(-3.0, 0.0);
  c.n0.resize(n);
  for (int i = 0; i < n; ++i)
    c.n0(i) = std::pow(10.0, logN0(rng));

  return c;
}

}  // namespace cram_test
