#pragma once
//
// Chebyshev Rational Approximation Method (CRAM) solver for the matrix
// exponential applied to a vector:  n(t) = exp(A t) n0.
//
// Implements the Incomplete Partial Factorization (IPF) form from
//   M. Pusa, "Higher-Order Chebyshev Rational Approximation Method and
//   Application to Burnup Equations," Nucl. Sci. Eng. 182:3, 297-318 (2016),
//   https://doi.org/10.13182/NSE15-26
//
//     exp(At) n0 ~ alpha0 * PROD_l ( I + 2 Re( alpha_l (At - theta_l I)^-1 ) ) n0
//
// applied incrementally to the running vector. Each pole requires one complex
// sparse linear solve; CRAM16 uses 8 poles, CRAM48 uses 24.
//
#include <Eigen/SparseCore>

namespace cram {

enum class CramOrder { CRAM16 = 16, CRAM48 = 48 };

// Solve n(t) = exp(A * dt) * n0.
// `A` is the real burnup matrix (see chain.hpp for the sign convention).
// CRAM48 is the safe default; CRAM16 is faster and usually accurate to ~1e-5
// relative on realistic depletion problems.
Eigen::VectorXd cramSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& n0,
                          double dt, CramOrder order = CramOrder::CRAM48);

}  // namespace cram
