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

void CramSolver::prepare(const Eigen::SparseMatrix<double>& A, double dt) {
  if (A.rows() != A.cols())
    throw std::invalid_argument("cram: CramSolver::prepare requires square A");

  n_ = static_cast<int>(A.rows());
  prepared_ = false;

  const Eigen::SparseMatrix<cd> Adt = (A.cast<cd>() * cd(dt, 0.0)).eval();
  Eigen::SparseMatrix<cd> I(n_, n_);
  I.setIdentity();

  // The sparsity of (A*dt - theta_l*I) is identical for every pole, but each
  // pole holds its own SparseLU instance because the numerical factorization
  // depends on theta_l. analyzePattern still costs negligibly compared to
  // factorize, and keeping per-pole solvers means apply() never has to
  // rebuild M.
  for (std::size_t l = 0; l < theta_.size(); ++l) {
    Eigen::SparseMatrix<cd> M = Adt - (theta_[l] * I);
    M.makeCompressed();
    auto& lu = *lus_[l];
    lu.analyzePattern(M);
    lu.factorize(M);
    if (lu.info() != Eigen::Success)
      throw std::runtime_error("cram: SparseLU factorization failed in CramSolver::prepare");
  }
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
