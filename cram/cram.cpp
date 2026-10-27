#include "cram/cram.hpp"

#include <Eigen/SparseLU>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "cram/cram_poles.hpp"
#include "cram/cram_solver.hpp"
#include "cram/fill_ordering_internal.hpp"

namespace cram {
namespace {

using cd = std::complex<double>;

Eigen::VectorXd ipfCram(const Eigen::SparseMatrix<double>& Araw, const Eigen::VectorXd& n0,
                        double dt, const CramPoles& poles) {
  const int N = static_cast<int>(n0.size());
  Eigen::SparseMatrix<cd> I(N, N);
  I.setIdentity();

  // Factorize P (A dt - theta I) P^T for the fill-reducing ordering P (see
  // fill_ordering_internal.hpp), carrying the vector into that ordering and
  // back: exp(A dt) = P^T exp(P A dt P^T) P.
  const Eigen::SparseMatrix<cd> scaled = (Araw.cast<cd>() * cd(dt, 0.0)).eval();
  Eigen::SparseMatrix<cd> M = scaled - (poles.theta[0] * I);
  M.makeCompressed();
  const internal::SymmetricPermutation P = internal::fillReducingOrder(M);
  const Eigen::SparseMatrix<cd> A = P * scaled * P.transpose();

  Eigen::VectorXd y = P * n0;  // stays real: we add 2*Re(...) each pole
  Eigen::SparseLU<Eigen::SparseMatrix<cd>, internal::PresetOrdering<int>> lu;

  // The sparsity pattern of (A - theta_l * I) does not depend on the pole, so
  // run the symbolic analysis once and only refactorize numerically per pole.
  M = A - (poles.theta[0] * I);
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
  return P.transpose() * (y * poles.alpha0);
}

}  // namespace

void requireTimeStep(double dt) {
  if (!std::isfinite(dt) || dt < 0.0)
    throw std::invalid_argument("cram: the time step must be finite and non-negative, got " +
                                std::to_string(dt));
}

Eigen::VectorXd cramSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& n0,
                          double dt, CramOrder order) {
  if (A.rows() != A.cols() || A.rows() != n0.size())
    throw std::invalid_argument("cram: A must be square with rows() == n0.size()");
  requireTimeStep(dt);

  // exp(A * 0) = I, and exp(0 * dt) = I — return n0 unchanged either way.
  if (dt == 0.0 || A.nonZeros() == 0)
    return n0;

  return ipfCram(A, n0, dt, cramPoles(order));
}

}  // namespace cram
