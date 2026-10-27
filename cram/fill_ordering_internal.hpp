#pragma once
//
// Internal fill-reducing ordering shared by cramSolve (one-shot) and CramSolver
// (cached). NOT part of the public API -- do not include from anywhere outside
// the cram library implementation.
//
#include <Eigen/SparseCore>

namespace cram::internal {

using SymmetricPermutation = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;

// A symmetric permutation P under which P*A*P^T is block lower triangular, with
// COLAMD ordering the rows and columns inside each diagonal block. Factorize
// the permuted matrix without a column ordering of its own: the ordering is already
// done, and a solver that reordered columns again would undo it -- use
// PresetOrdering (cram_solver.hpp) for that.
//
// This is the structure KLU exploits. A burnup matrix is far from irreducible
// -- decay and capture mostly run one way through the chain, and the cycles
// come only from (n,2n)-style loops and fission feeding back -- so its strongly
// connected components are small, the blocks below the diagonal produce no
// fill at all, and COLAMD is left a much smaller problem in each block. A
// decay-only matrix is acyclic, every block is a single nuclide, and the
// factorization has no fill whatever.
//
// Measured against SparseLU's own COLAMD over the whole matrix, at N = 1675
// with CRAM48 (MSVC, Release). Burnup benchmark shape: L+U shrinks from 2.24x
// to 1.88x nnz(A - theta*I), CramSolver::prepare() from 31 to 16 ms, apply()
// grows from 0.71 to 0.78 ms, and a 100-step predictor march drops 26%. Acyclic
// decay shape: 1.70x -> 1.40x (no fill-in at all), prepare() 19.5 -> 9.3 ms,
// apply() 0.62 -> 0.55 ms. The OpenMC simplified PWR chain (N = 228) changes
// little in fill (1.15x -> 1.13x) but factorizes ~20% faster.
//
// Depends only on the sparsity pattern of A (entries present, whatever their
// value), and is deterministic: a function of that pattern alone.
template <typename Scalar>
SymmetricPermutation fillReducingOrder(const Eigen::SparseMatrix<Scalar>& A);

}  // namespace cram::internal
