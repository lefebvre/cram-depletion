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

void CramSolver::prepare(const Eigen::SparseMatrix<double>& A, double dt) {
  if (A.rows() != A.cols())
    throw std::invalid_argument("cram: CramSolver::prepare requires square A");

  n_ = static_cast<int>(A.rows());
  prepared_ = false;

  const Eigen::SparseMatrix<cd> Adt = (A.cast<cd>() * cd(dt, 0.0)).eval();
  Eigen::SparseMatrix<cd> I(n_, n_);
  I.setIdentity();

  // The sparsity of (A*dt - theta_l*I) is identical for every pole, so one
  // pole matrix settles the pattern question for all of them. Compare against
  // the last prepared pattern before touching any solver: if it still holds,
  // every pole's symbolic analysis is still valid and only the numeric
  // factorization has to be redone. On a many-region sweep, where each region
  // shares the chain's topology and differs only in reaction rates, that is
  // the difference between analyzing once and analyzing once per region.
  Eigen::SparseMatrix<cd> M0 = Adt - (theta_[0] * I);
  M0.makeCompressed();
  const bool reuse = symbolicValid_ && patternMatches(M0);

  // A failed factorization leaves the symbolic state untrustworthy, so drop it
  // up front and only reinstate it once every pole has gone through cleanly.
  symbolicValid_ = false;

  for (std::size_t l = 0; l < theta_.size(); ++l) {
    // Each pole keeps its own SparseLU because the numeric factorization
    // depends on theta_l; only the symbolic half is shared.
    Eigen::SparseMatrix<cd> M = (l == 0) ? M0 : (Adt - (theta_[l] * I));
    if (l != 0)
      M.makeCompressed();
    auto& lu = *lus_[l];
    if (!reuse)
      lu.analyzePattern(M);
    lu.factorize(M);
    if (lu.info() != Eigen::Success)
      throw std::runtime_error("cram: SparseLU factorization failed in CramSolver::prepare");
  }

  if (!reuse)
    recordPattern(M0);
  symbolicValid_ = true;
  reusedSymbolic_ = reuse;
  prepared_ = true;
}

Eigen::VectorXd CramSolver::apply(const Eigen::VectorXd& n0) const {
  if (!prepared_)
    throw std::logic_error("cram: CramSolver::apply called before prepare");
  if (n0.size() != n_)
    throw std::invalid_argument("cram: CramSolver::apply n0.size() does not match prepared A");

  Eigen::VectorXd y = n0;
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
