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
  // Throws std::invalid_argument if A is not square; std::runtime_error if
  // any pole's factorization fails (singular matrix, NaNs in A, etc.).
  void prepare(const Eigen::SparseMatrix<double>& A, double dt);

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

  CramOrder order_;
  int n_ = 0;
  double alpha0_ = 0.0;
  std::vector<cd> theta_;
  std::vector<cd> alpha_;
  // One SparseLU per pole. We can't keep these in a std::vector of value
  // type because Eigen::SparseLU is non-copyable / non-movable in the
  // assignment sense we'd want — use unique_ptr to give us stable storage.
  std::vector<std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<cd>>>> lus_;
  bool prepared_ = false;
};

}  // namespace cram
