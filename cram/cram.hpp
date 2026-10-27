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
// sparse linear solve; CRAM16 uses 8 poles, CRAM48 uses 24. The coefficients,
// and CramOrder itself, are published in cram/cram_poles.hpp.
//
#include <Eigen/SparseCore>

#include "cram/cram_poles.hpp"

namespace cram {

// Solve n(t) = exp(A * dt) * n0.
// `A` is the real burnup matrix (see chain.hpp for the sign convention).
// CRAM48 is the safe default; CRAM16 is faster and usually accurate to ~1e-5
// relative on realistic depletion problems. Throws std::invalid_argument for a
// non-square A, a size mismatch, or a dt rejected by requireTimeStep().
Eigen::VectorXd cramSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& n0,
                          double dt, CramOrder order = CramOrder::CRAM48);

// Throws std::invalid_argument unless dt is finite and non-negative. CRAM's
// rational approximation is accurate only for eigenvalues of A*dt near the
// negative real axis; a burnup matrix's eigenvalues have non-positive real
// part, so a negative dt reflects them onto the positive side, where the
// approximation returns a wrong answer rather than an error. Every entry point
// that takes a time step applies this check.
void requireTimeStep(double dt);

}  // namespace cram
