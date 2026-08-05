#pragma once
//
// Stateful CRAM solver that caches the per-pole LU factorizations of
// (A*dt - theta_l*I). Use when the same (A, dt) pair is reused across many
// matrix-exponential evaluations — decay-only marches over many sub-steps,
// or intra-step sub-stepping inside one burnup step.
//
// For a one-shot evaluation, use the free function cramSolve() in cram.hpp.
//
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <complex>
#include <memory>
#include <vector>

#include "cram/cram.hpp"

namespace cram {

class CramSolver {
public:
  explicit CramSolver(CramOrder order = CramOrder::CRAM48);

  // Factorize (A*dt - theta_l*I) for every pole l. Must be called at least
  // once before apply(). Calling prepare() again replaces the cache —
  // there is no auto-detection of "same arguments," so the caller is
  // responsible for only re-preparing when A or dt actually changes.
  //
  // The symbolic analysis IS reused automatically. Each pole's factorization
  // splits into a symbolic phase that depends only on sparsity and a numeric
  // phase that depends on the values; the symbolic phase is ~40% of prepare()
  // at N~1700 and is identical for every pole and for every matrix with the
  // same pattern. When a later prepare() is handed a matrix whose pattern
  // matches the previous one, that phase is skipped.
  //
  // This is what makes the many-region workload affordable: a reactor model
  // with one burnup matrix per depletion region has the same chain topology
  // — and therefore the same sparsity — in every region, so the symbolic
  // analysis is genuinely needed once, not once per region. Hold one
  // CramSolver per worker thread and call prepare() on it for each region in
  // turn; only the numeric factorization is repaid. Pattern equality is
  // checked by comparing the compressed index arrays, which is O(nnz) and far
  // cheaper than the analysis it avoids. A mismatch is not an error: it simply
  // falls back to a full analysis.
  //
  // Throws std::invalid_argument if A is not square; std::runtime_error if
  // any pole's factorization fails (singular matrix, NaNs in A, etc.).
  void prepare(const Eigen::SparseMatrix<double>& A, double dt);

  // True if the most recent prepare() was able to reuse the symbolic analysis
  // from the one before it. Exposed for tests and for callers who want to
  // confirm their region loop is actually hitting the fast path.
  bool reusedSymbolicAnalysis() const { return reusedSymbolic_; }

  // Apply the cached exp(A*dt) operator to n0. n0.size() must match the
  // size of A passed to the most recent prepare().
  //
  // Throws std::logic_error if prepare() has not been called;
  // std::invalid_argument if n0.size() is wrong; std::runtime_error if a
  // per-pole solve reports failure.
  Eigen::VectorXd apply(const Eigen::VectorXd& n0) const;

  int size() const { return n_; }
  CramOrder order() const { return order_; }
  bool prepared() const { return prepared_; }

  // -- Thread-safety ------------------------------------------------------
  //
  // A single CramSolver instance is NOT thread-safe. prepare() mutates the
  // cache; apply() reads cached Eigen::SparseLU state that is not safe to
  // share across threads (the underlying SuperNodalMatrix workspace is
  // touched on every solve). To parallelize, construct one CramSolver per
  // thread — the pole tables themselves are immutable shared constants, so
  // there is no contention between instances.

private:
  using cd = std::complex<double>;
  using StorageIndex = Eigen::SparseMatrix<cd>::StorageIndex;

  // Sparsity of the pole matrices from the last successful prepare(), kept so
  // the next one can tell whether the symbolic analysis still applies. Cheap
  // to hold: O(n + nnz) indices against K stored factorizations.
  bool patternMatches(const Eigen::SparseMatrix<cd>& m) const;
  void recordPattern(const Eigen::SparseMatrix<cd>& m);

  CramOrder order_;
  int n_ = 0;
  double alpha0_ = 0.0;
  std::vector<cd> theta_;
  std::vector<cd> alpha_;
  std::vector<StorageIndex> patternOuter_;
  std::vector<StorageIndex> patternInner_;
  bool symbolicValid_ = false;
  bool reusedSymbolic_ = false;
  // One SparseLU per pole. We can't keep these in a std::vector of value
  // type because Eigen::SparseLU is non-copyable / non-movable in the
  // assignment sense we'd want — use unique_ptr to give us stable storage.
  std::vector<std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<cd>>>> lus_;
  bool prepared_ = false;
};

}  // namespace cram
