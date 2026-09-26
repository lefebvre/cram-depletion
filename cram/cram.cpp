#include "cram/cram.hpp"

#include <Eigen/SparseLU>
#include <complex>
#include <cstddef>
#include <stdexcept>

#include "cram/cram_poles.hpp"

namespace cram {
namespace {

using cd = std::complex<double>;

Eigen::VectorXd ipfCram(const Eigen::SparseMatrix<double>& Araw, const Eigen::VectorXd& n0,
                        double dt, const CramPoles& poles) {
  const int N = static_cast<int>(n0.size());
  const Eigen::SparseMatrix<cd> A = (Araw.cast<cd>() * cd(dt, 0.0)).eval();
  Eigen::SparseMatrix<cd> I(N, N);
  I.setIdentity();

  Eigen::VectorXd y = n0;  // stays real: we add 2*Re(...) each pole
  Eigen::SparseLU<Eigen::SparseMatrix<cd>> lu;

  // The sparsity pattern of (A - theta_l * I) does not depend on the pole, so
  // run the symbolic analysis once and only refactorize numerically per pole.
  Eigen::SparseMatrix<cd> M = A - (poles.theta[0] * I);
  M.makeCompressed();
  lu.analyzePattern(M);

  for (std::size_t l = 0; l < poles.theta.size(); ++l) {
    if (l > 0) {
      M = A - (poles.theta[l] * I);
      M.makeCompressed();
    }
    lu.factorize(M);
    if (lu.info() != Eigen::Success)
      throw std::runtime_error("cram: SparseLU factorization failed");
    Eigen::VectorXcd x = lu.solve(y.cast<cd>());
    if (lu.info() != Eigen::Success)
      throw std::runtime_error("cram: SparseLU solve failed");
    y += 2.0 * (poles.alpha[l] * x).real();
  }
  return y * poles.alpha0;
}

}  // namespace

Eigen::VectorXd cramSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& n0,
                          double dt, CramOrder order) {
  if (A.rows() != A.cols() || A.rows() != n0.size())
    throw std::invalid_argument("cram: A must be square with rows() == n0.size()");

  // exp(A * 0) = I, and exp(0 * dt) = I — return n0 unchanged either way.
  if (dt == 0.0 || A.nonZeros() == 0)
    return n0;

  return ipfCram(A, n0, dt, cramPoles(order));
}

}  // namespace cram
