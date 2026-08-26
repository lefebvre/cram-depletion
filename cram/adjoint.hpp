#pragma once
//
// Adjoint depletion for linear systems, and the piecewise-constant marches it
// pairs with.
//
// For dn/dt = A n with A independent of n, a response R = <w, n(T)> on the
// final composition has the adjoint (importance) n*(t) satisfying
//
//     dn*/dt = -A^T n*,   n*(T) = w,
//
// so that <n*(t), n(t)> is constant in t and n*(0) = dR/dn0: the sensitivity
// of the response to every initial number density at once, from one backward
// solve. Because exp(A^T dt) = exp(A dt)^T, the backward step is the same CRAM
// solve on the transposed matrix, and <w, exp(A dt) n0> = <exp(A^T dt) w, n0>
// is the duality identity every result here must satisfy.
//
// Scope: A must be composition-independent on each interval -- decay-only,
// constant flux, or a piecewise-constant flux schedule (one matrix per
// interval). Constant-power normalization makes A depend on n through the
// flux, and its adjoint needs the rank-one feedback Jacobian and a discrete
// adjoint of the chosen integrator; that is not implemented, and
// intervalMatrices() refuses a system in that mode rather than freezing the
// flux silently.
//
#include <Eigen/SparseCore>
#include <vector>

#include "cram/cram.hpp"
#include "cram/deplete.hpp"
#include "cram/integrator.hpp"

namespace cram {

// A^T as a compressed matrix. Eigen's A.transpose() is an expression; the
// solvers want a materialized SparseMatrix.
Eigen::SparseMatrix<double> transposed(const Eigen::SparseMatrix<double>& A);

// n*(t - dt) = exp(A^T dt) w: one backward adjoint step. Equal to
// cramSolve(transposed(A), w, dt, order).
Eigen::VectorXd cramSolveAdjoint(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& w,
                                 double dt, CramOrder order = CramOrder::CRAM48);

// Forward march through per-interval matrices: n_{k+1} = exp(A_k dt_k) n_k.
// A.size() == dts.size(); the result has dts.size() + 1 points. Every interval
// is one CramSolver::prepare() (the symbolic analysis is reused across
// intervals, since they share a pattern) and one apply(). This is the
// linear-system counterpart of deplete(); the two agree whenever the
// integrator's builder returns A_k on interval k.
DepletionResult depleteLinear(const std::vector<Eigen::SparseMatrix<double>>& A,
                              const std::vector<double>& dts, const Eigen::VectorXd& n0,
                              CramOrder order = CramOrder::CRAM48);

// Adjoint trajectory, stored forward in time so nStar[k] sits at the same
// point as DepletionResult::n[k]: nStar[K] = w and nStar[k] = exp(A_k^T dt_k)
// nStar[k+1]. nStar.front() is dR/dn0 for R = <w, n(T)>.
struct AdjointResult {
  std::vector<double> time;
  std::vector<Eigen::VectorXd> nStar;
};

AdjointResult adjointDeplete(const std::vector<Eigen::SparseMatrix<double>>& A,
                             const std::vector<double>& dts, const Eigen::VectorXd& w,
                             CramOrder order = CramOrder::CRAM48);

// One burnup matrix per interval for a fixed-flux system: interval k is
// assembled at flux[k]. The system must be in constant-flux mode (its flux is
// overwritten by each interval's value and left at the last one); a
// constant-power system throws std::invalid_argument, since its matrix depends
// on the composition and no single A_k represents the interval.
std::vector<Eigen::SparseMatrix<double>> intervalMatrices(DepletionSystem& system,
                                                          const std::vector<double>& flux);

}  // namespace cram
