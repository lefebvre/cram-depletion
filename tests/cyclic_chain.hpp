#pragma once
//
// A seeded, deterministic synthetic burnup matrix that is *cyclic*, unlike the
// linear chains used by the analytic-Bateman tests.
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
#include <cmath>
#include <random>
#include <vector>

namespace cram_test {

struct CyclicChain {
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd n0;
};

// Build an n x n cyclic burnup matrix. `seed` fixes the topology and the rates,
// so the same (n, seed) always yields a bit-identical matrix on every platform:
// std::mt19937 and the integer/real distributions used here are specified by the
// standard, and no floating-point value depends on iteration order.
inline CyclicChain buildCyclicChain(int n, unsigned seed = 20260804u) {
  // Seeded through std::seed_seq rather than mt19937's single-value constructor.
  // That constructor fills all 624 state words from one 32-bit value, so only
  // 2^32 of the engine's 2^19937 initial states are reachable and the resulting
  // state is structured rather than well mixed. seed_seq's expansion reaches
  // ~2^256 states here and mixes them, at no runtime cost.
  //
  // Reproducibility is unaffected, which is what this fixture actually requires:
  // both seed_seq's generate() and mt19937 are fully specified by the standard,
  // so a given `seed` yields a bit-identical matrix on every implementation --
  // which is what lets the frozen vectors in cram_golden_data.hpp be portable.
  //
  // The fixed words below are arbitrary mixing constants (golden ratio and
  // xxHash primes); they pad the entropy without making the stream depend on
  // anything but `seed`. Changing them, or the seeding procedure, changes every
  // generated matrix and requires regenerating cram_golden_data.hpp.
  std::seed_seq sequence{seed,        0x9e3779b9u, 0x85ebca6bu, 0xc2b2ae35u,
                         0x27d4eb2fu, 0x165667b1u, 0xd3a2646cu, 0xfd7046c5u};
  std::mt19937 rng(sequence);
  std::uniform_real_distribution<double> logLambda(-9.0, -1.0);  // 1e-9 .. 1e-1 /s
  std::uniform_int_distribution<int> nModes(1, 3);
  std::uniform_int_distribution<int> hop(1, 4);
  std::uniform_real_distribution<double> branch(0.15, 1.0);
  std::uniform_int_distribution<int> stableRoll(0, 19);  // ~5% stable
  std::uniform_int_distribution<int> lowIndex(0, (n / 4) - 1);

  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n) * 5);

  for (int i = 0; i < n; ++i) {
    // ~5% of nuclides are stable sinks: no diagonal, no production. Index 0 is
    // always active so the seeded inventory actually moves.
    if (i != 0 && stableRoll(rng) == 0)
      continue;

    const double lambda = std::pow(10.0, logLambda(rng));
    // std::size_t because `modes` and `k` only ever count and index the two
    // std::vectors below. The nuclide indices (n, i, j) stay signed: `j` is
    // deliberately negative in flight on the up-chain branch, `n` is compared
    // against it, and Eigen's Index and Triplet StorageIndex are signed too.
    const std::size_t modes = static_cast<std::size_t>(nModes(rng));

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
      if (j == i)  // never fold a branch into the diagonal
        j = (i + 1 < n) ? i + 1 : i - 1;
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
