#include "cram/cram_solver.hpp"

#include <algorithm>
#include <stdexcept>

#include "cram/cram_poles_internal.hpp"

namespace cram {

CramSolver::CramSolver(CramOrder order) : order_(order) {
  using namespace cram::internal;
  if (order == CramOrder::CRAM16) {
    alpha0_ = kAlpha0_16;
    theta_.assign(kTheta16.begin(), kTheta16.end());
    alpha_.assign(kAlpha16.begin(), kAlpha16.end());
  } else {
    alpha0_ = kAlpha0_48;
    theta_.assign(kTheta48.begin(), kTheta48.end());
    alpha_.assign(kAlpha48.begin(), kAlpha48.end());
  }
  lus_.resize(theta_.size());
  std::generate(lus_.begin(), lus_.end(),
                [] { return std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<cd>>>(); });
}

bool CramSolver::patternMatches(const Eigen::SparseMatrix<cd>& m) const {
  if (patternOuter_.empty())
    return false;
  if (static_cast<std::size_t>(m.outerSize()) + 1 != patternOuter_.size())
    return false;
  if (static_cast<std::size_t>(m.nonZeros()) != patternInner_.size())
    return false;
  return std::equal(patternOuter_.begin(), patternOuter_.end(), m.outerIndexPtr()) &&
         std::equal(patternInner_.begin(), patternInner_.end(), m.innerIndexPtr());
}

void CramSolver::recordPattern(const Eigen::SparseMatrix<cd>& m) {
  patternOuter_.assign(m.outerIndexPtr(), m.outerIndexPtr() + m.outerSize() + 1);
  patternInner_.assign(m.innerIndexPtr(), m.innerIndexPtr() + m.nonZeros());
}

void CramSolver::beginPrepare(const Eigen::SparseMatrix<double>& A, double dt) {
  if (A.rows() != A.cols())
    throw std::invalid_argument("cram: CramSolver::prepare requires square A");

  n_ = static_cast<int>(A.rows());
  prepared_ = false;

  scaled_ = (A.cast<cd>() * cd(dt, 0.0)).eval();

  // The sparsity of (A*dt - theta_l*I) is identical for every pole, so one
  // pole matrix settles the pattern question for all of them. Decide reuse
  // here, before any pole is touched, so the pole loop reads a value nobody
  // writes -- which is what lets those calls run concurrently.
  Eigen::SparseMatrix<cd> I(n_, n_);
  I.setIdentity();
  Eigen::SparseMatrix<cd> M0 = scaled_ - (theta_[0] * I);
  M0.makeCompressed();
  reuseSymbolic_ = symbolicValid_ && patternMatches(M0);
  if (!reuseSymbolic_)
    recordPattern(M0);

  // A failed factorization leaves the symbolic state untrustworthy, so drop it
  // up front and only reinstate it in endPrepare() once every pole is clean.
  symbolicValid_ = false;
  poleOk_.assign(theta_.size(), 0);
  preparing_ = true;
}

void CramSolver::preparePole(std::size_t pole) {
  if (!preparing_)
    throw std::logic_error("cram: CramSolver::preparePole called outside begin/endPrepare");
  if (pole >= theta_.size())
    throw std::out_of_range("cram: CramSolver::preparePole index out of range");

  // Built locally so concurrent calls share nothing writable. Each pole keeps
  // its own SparseLU because the numeric factorization depends on theta_l;
  // only the symbolic half is common to all of them.
  Eigen::SparseMatrix<cd> I(n_, n_);
  I.setIdentity();
  Eigen::SparseMatrix<cd> M = scaled_ - (theta_[pole] * I);
  M.makeCompressed();

  auto& lu = *lus_[pole];
  if (!reuseSymbolic_)
    lu.analyzePattern(M);
  lu.factorize(M);
  // Reported rather than thrown: an exception escaping a caller's parallel
  // region would terminate. endPrepare() raises it on the calling thread.
  poleOk_[pole] = (lu.info() == Eigen::Success) ? 1u : 0u;
}

void CramSolver::endPrepare() {
  if (!preparing_)
    throw std::logic_error("cram: CramSolver::endPrepare called without beginPrepare");
  preparing_ = false;

  for (std::size_t l = 0; l < poleOk_.size(); ++l) {
    if (poleOk_[l] == 0u) {
      scaled_.resize(0, 0);
      throw std::runtime_error("cram: SparseLU factorization failed in CramSolver::prepare");
    }
  }

  // The scaled matrix is only needed while poles are being factorized. Dropping
  // it here keeps a per-thread solver's footprint to the K factorizations.
  scaled_.resize(0, 0);
  symbolicValid_ = true;
  prepared_ = true;
}

void CramSolver::prepare(const Eigen::SparseMatrix<double>& A, double dt) {
  beginPrepare(A, dt);
  for (std::size_t l = 0; l < theta_.size(); ++l)
    preparePole(l);
  endPrepare();
}

Eigen::VectorXd CramSolver::apply(const Eigen::VectorXd& n0) const {
  if (!prepared_)
    throw std::logic_error("cram: CramSolver::apply called before prepare");
  if (n0.size() != n_)
    throw std::invalid_argument("cram: CramSolver::apply n0.size() does not match prepared A");

  Eigen::VectorXd y = n0;
  // This loop allocates two n-vectors per pole. Hoisting them into reusable
  // member scratch was tried and measured no faster at N = 256, 1024 or 1675 --
  // the triangular solves dominate so thoroughly that the allocations do not
  // show above run-to-run noise -- so the simpler code stands.
  for (std::size_t l = 0; l < theta_.size(); ++l) {
    const auto& lu = *lus_[l];
    Eigen::VectorXcd x = lu.solve(y.cast<cd>());
    if (lu.info() != Eigen::Success)
      throw std::runtime_error("cram: SparseLU solve failed in CramSolver::apply");
    y += 2.0 * (alpha_[l] * x).real();
  }
  return y * alpha0_;
}

}  // namespace cram
