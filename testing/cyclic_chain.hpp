#pragma once
//
// A seeded, deterministic synthetic burnup matrix that is *cyclic*, unlike the
// linear chains used by the analytic-Bateman tests.
//
// Shared by the golden regression test and the benchmark suite on purpose: what
// the benchmarks time and what the golden vector locks should be the same matrix
// structure, so the generator lives here rather than being duplicated. Changing
// it changes the frozen data in tests/cram_golden_data.hpp -- the golden test
// will fail if the topology or the RNG draw sequence shifts, which is intended.
//
// Why this exists: every other solver test builds a strictly lower-triangular A
// (a linear decay chain), for which Eigen::SparseLU does no pivoting and
// produces no fill-in. That leaves the pole-matrix assembly path exercised only
// in its easiest case. Real burnup matrices have cycles -- (n,gamma) walks the
// chain up in A, decay walks it back down, and fission products link high
// indices back to low ones -- so the factorization has substantial fill-in.
// This generator reproduces that structure.
//
// Two properties are guaranteed by construction and are relied on by the tests:
//
//   1. Every column sums to zero. Each nuclide's branching fractions sum to
//      exactly 1.0 (the last branch is computed as 1 - sum(others)), and every
//      daughter lands inside [0, n), so no atoms leave the tracked set.
//      Total atom count is therefore invariant under exp(A*dt) to within
//      rounding, at any dt.
//
//   2. A fraction of the nuclides are stable: they emit no triplets at all, so
//      their diagonal entry is structurally absent from A and (A*dt - theta*I)
//      must create it. That mirrors the terminal-stable-nuclide case the
//      linear-chain tests cover, but at N where fill-in is also in play.
//
#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace cram_test {

struct CyclicChain {
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd n0;
};

// Build an n x n cyclic burnup matrix. `seed` fixes the topology and the rates,
// so the same (n, seed) yields a bit-identical matrix on every toolchain CI
// exercises, which is what lets the frozen vectors in cram_golden_data.hpp be
// checked in rather than regenerated per platform. No floating-point value here
// depends on iteration order.
//
// That portability is verified rather than guaranteed by the standard; see the
// note on the distributions below for what is and is not promised.
inline CyclicChain buildCyclicChain(int n, unsigned seed = 20260804u) {
  // Seeded through std::seed_seq rather than mt19937's single-value constructor.
  // That constructor fills all 624 state words from one 32-bit value, so only
  // 2^32 of the engine's 2^19937 initial states are reachable and the resulting
  // state is structured rather than well mixed. seed_seq's expansion reaches
  // ~2^256 states here and mixes them, at no runtime cost.
  //
  // Reproducibility is unaffected, which is what this fixture actually requires:
  // seed_seq::generate() and mt19937 are both algorithm-specified, so the engine
  // emits the same word sequence from a given `seed` on every implementation.
  //
  // The distributions below carry no such promise. [rand.dist.uni.int] and
  // [rand.dist.uni.real] specify the distribution a generator produces, not the
  // algorithm that produces it, so an implementation may consume a different
  // number of engine words or map them differently -- and libstdc++, libc++ and
  // the MSVC STL are known to differ in practice. A fixed engine sequence is
  // therefore necessary for a portable golden vector but not sufficient.
  //
  // What actually makes cram_golden_data.hpp portable is that CI checks it: the
  // golden test runs under libstdc++ and the MSVC STL on every push, so a
  // divergence shows up as a test failure rather than silently. libc++ is not in
  // the matrix; adding a macOS or -stdlib=libc++ job may well require
  // regenerating the vectors, and the golden test is what will say so.
  //
  // The fixed words below are arbitrary mixing constants (golden ratio and
  // xxHash primes); they pad the entropy without making the stream depend on
  // anything but `seed`. Changing them, or the seeding procedure, changes every
  // generated matrix and requires regenerating cram_golden_data.hpp.
  std::seed_seq sequence{seed,        0x9e3779b9u, 0x85ebca6bu, 0xc2b2ae35u,
                         0x27d4eb2fu, 0x165667b1u, 0xd3a2646cu, 0xfd7046c5u};
  std::mt19937 rng(sequence);
  std::uniform_real_distribution<double> logLambda(-9.0, -1.0);  // 1e-9 .. 1e-1 /s
  // std::size_t rather than the default int: this draw feeds `modes` below,
  // which only ever counts and indexes std::vectors. The result type does not
  // perturb the stream -- for a range this narrow libstdc++ and the MSVC STL
  // both consume one engine word per attempt and reject on the same threshold
  // whatever its width -- so the golden vectors are unaffected. That is
  // measured, not promised; see the note above on what the standard specifies.
  std::uniform_int_distribution<std::size_t> nModes(1, 3);
  std::uniform_int_distribution<int> hop(1, 4);
  std::uniform_real_distribution<double> branch(0.15, 1.0);
  std::uniform_int_distribution<int> stableRoll(0, 19);  // ~5% stable
  // The upper bound is clamped so the range stays valid at n < 4, where
  // n/4 - 1 is -1: uniform_int_distribution requires a <= b, so that would be
  // UB here at construction -- before any nuclide is drawn, and therefore even
  // for n == 0, whose loops never run. For n >= 4 the range is exactly
  // n/4 - 1 as before, so the draw sequence and the golden data are untouched.
  std::uniform_int_distribution<int> lowIndex(0, std::max(1, n / 4) - 1);

  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n) * 5);

  for (int i = 0; i < n; ++i) {
    // ~5% of nuclides are stable sinks: no diagonal, no production. Index 0 is
    // always active so the seeded inventory actually moves.
    if (i != 0 && stableRoll(rng) == 0)
      continue;

    const double lambda = std::pow(10.0, logLambda(rng));
    // `modes` and `k` only ever count and index the two std::vectors below, so
    // both are std::size_t. The nuclide indices (n, i, j) stay signed: `j` is
    // deliberately negative in flight on the up-chain branch, `n` is compared
    // against it, and Eigen's Index and Triplet StorageIndex are signed too.
    const std::size_t modes = nModes(rng);

    // Daughters: alternate down-chain (decay) and up-chain ((n,gamma)-like) so
    // the matrix has entries on both sides of the diagonal, and every 16th
    // nuclide also gets a long-range link back to a low index (fission-like).
    // Long-range back-edges are what force real fill-in.
    std::vector<int> targets;
    targets.reserve(modes);
    for (std::size_t k = 0; k < modes; ++k) {
      int j;
      if (i % 8 == 7 && k == 0) {
        j = lowIndex(rng);  // long-range back-edge
      } else if (k % 2 == 0) {
        j = i + hop(rng);  // down-chain
      } else {
        j = i - hop(rng);  // up-chain
      }
      if (j < 0)
        j = 0;
      if (j >= n)
        j = n - 1;
      // Never fold a branch into the diagonal: step up if there is room, else
      // down. At n == 1 there is nowhere to go and the branch has to land on
      // the diagonal -- taking i - 1 there would emit row -1. The column still
      // sums to zero there, since the branch cancels the removal term.
      if (j == i)
        j = (i + 1 < n) ? i + 1 : (i > 0 ? i - 1 : i);
      targets.push_back(j);
    }

    // Branching fractions summing to exactly 1.0: draw all but the last, then
    // take the remainder. Guards against a negative remainder by falling back
    // to an even split.
    // `modes` is drawn from [1, 3], so `modes - 1` cannot wrap.
    std::vector<double> br(modes);
    const double modesAsDouble = static_cast<double>(modes);
    double drawn = 0.0;
    for (std::size_t k = 0; k + 1 < modes; ++k) {
      br[k] = branch(rng) / modesAsDouble;
      drawn += br[k];
    }
    if (drawn >= 1.0) {
      for (std::size_t k = 0; k < modes; ++k)
        br[k] = 1.0 / modesAsDouble;
    } else {
      br[modes - 1] = 1.0 - drawn;
    }

    t.emplace_back(i, i, -lambda);  // total removal
    for (std::size_t k = 0; k < modes; ++k)
      t.emplace_back(targets[k], i, lambda * br[k]);
  }

  CyclicChain c;
  c.A.resize(n, n);
  c.A.setFromTriplets(t.begin(), t.end());  // sums duplicates
  c.A.makeCompressed();

  // Seed every nuclide, spread over three decades. A sparse initial inventory
  // leaves most of the solution vector at exactly zero at any dt -- unreachable
  // or underflowed -- which would make a golden vector a weak lock, since half
  // its entries would carry no information. A full inventory (which is also what
  // a real depletion problem looks like) keeps every component meaningful.
  std::uniform_real_distribution<double> logN0(-3.0, 0.0);
  c.n0.resize(n);
  for (int i = 0; i < n; ++i)
    c.n0(i) = std::pow(10.0, logN0(rng));

  return c;
}

}  // namespace cram_test
