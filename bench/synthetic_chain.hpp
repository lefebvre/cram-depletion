#pragma once
//
// Synthetic depletion chain used for benchmarking. Builds an N-nuclide
// burnup matrix whose sparsity profile is representative of a real depletion
// chain (each nuclide has 1-3 outgoing transitions to nearby indices, with
// branching ratios summing to one). Deterministic and reproducible.
//
// Why synthetic: loading a real ENDF/B chain at bench time would force every
// run to depend on a multi-MB data file. The numeric pattern doesn't need
// to be physical for performance measurement — only the shape of A does.
//
#include <Eigen/SparseCore>

namespace cram_bench {

struct SyntheticChain {
  Eigen::SparseMatrix<double> A;  // N x N burnup matrix, OpenMC sign convention
  Eigen::VectorXd n0;             // initial inventory (single seed at index 0)
};

// Build a chain with `n` nuclides. Decay constants span roughly eight orders
// of magnitude — the same kind of stiffness CRAM is designed to handle.
SyntheticChain buildSyntheticChain(int n, unsigned seed = 0xC8A48u);

}  // namespace cram_bench
