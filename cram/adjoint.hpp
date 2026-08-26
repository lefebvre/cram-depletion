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
// The same trajectories give first-order sensitivities to the matrix itself:
// dR/dA_ij = integral_0^T n*_i(t) n_j(t) dt, from which the sensitivity to any
// decay constant or cross section follows by the chain rule (see
// rateSensitivities() and the two contractions below).
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
#include <optional>
#include <vector>

#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "cram/reaction.hpp"

namespace cram {

// A^T as a compressed matrix. Eigen's A.transpose() is an expression; the
// solvers want a materialized SparseMatrix.
Eigen::SparseMatrix<double> transposed(const Eigen::SparseMatrix<double>& A);

// n*(t - dt) = exp(A^T dt) w: one backward adjoint step. Equal to
// cramSolve(transposed(A), w, dt, order).
Eigen::VectorXd cramSolveAdjoint(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& w,
                                 double dt, CramOrder order = CramOrder::CRAM48);

// --- Generators for inhomogeneous and integrated problems -------------------

// The (N+1) x (N+1) generator [[A, s], [0, 0]] of dn/dt = A n + s: applied to
// (n0; 1), its exponential carries n(t) in the first N entries and leaves the
// trailing 1 untouched. This is how a constant feed or removal source, or a
// time-integrated adjoint response, enters a plain matrix exponential.
Eigen::SparseMatrix<double> sourceGenerator(const Eigen::SparseMatrix<double>& A,
                                            const Eigen::VectorXd& s);

// The 2N x 2N generator [[A, 0], [I, 0]]: applied to (n0; 0), its exponential
// carries n(t) in the top block and the exact integral_0^t n(s) ds in the
// bottom one (Van Loan, 1978).
Eigen::SparseMatrix<double> augmentedGenerator(const Eigen::SparseMatrix<double>& A);

// --- Piecewise-constant marches -----------------------------------------------

// Forward march through per-interval matrices: n_{k+1} = exp(A_k dt_k) n_k.
// A.size() == dts.size(); the result has dts.size() + 1 points. Every interval
// is one CramSolver::prepare() (the symbolic analysis is reused across
// intervals, since they share a pattern) and one apply(). This is the
// linear-system counterpart of deplete(); the two agree whenever the
// integrator's builder returns A_k on interval k.
DepletionResult depleteLinear(const std::vector<Eigen::SparseMatrix<double>>& A,
                              const std::vector<double>& dts, const Eigen::VectorXd& n0,
                              CramOrder order = CramOrder::CRAM48);

// integral_0^T n(t) dt over the schedule, exact per interval through
// augmentedGenerator().
Eigen::VectorXd integratedInventory(const std::vector<Eigen::SparseMatrix<double>>& A,
                                    const std::vector<double>& dts, const Eigen::VectorXd& n0,
                                    CramOrder order = CramOrder::CRAM48);

// Adjoint trajectory, stored forward in time so nStar[k] sits at the same
// point as DepletionResult::n[k]. nStar.front() is dR/dn0.
//
// For a final-time response R = <w, n(T)>: nStar[K] = w and
// nStar[k] = exp(A_k^T dt_k) nStar[k+1]; integratedWeight is empty.
//
// For a time-integrated response R = integral_0^T <w, n(t)> dt: nStar[K] = 0,
// each backward step also accumulates the source w, and integratedWeight
// holds w so rateSensitivities() can evaluate the adjoint between points.
struct AdjointResult {
  std::vector<double> time;
  std::vector<Eigen::VectorXd> nStar;
  std::optional<Eigen::VectorXd> integratedWeight;
};

AdjointResult adjointDeplete(const std::vector<Eigen::SparseMatrix<double>>& A,
                             const std::vector<double>& dts, const Eigen::VectorXd& w,
                             CramOrder order = CramOrder::CRAM48);

AdjointResult adjointDepleteIntegrated(const std::vector<Eigen::SparseMatrix<double>>& A,
                                       const std::vector<double>& dts, const Eigen::VectorXd& w,
                                       CramOrder order = CramOrder::CRAM48);

// One burnup matrix per interval for a fixed-flux system: interval k is
// assembled at flux[k]. The system must be in constant-flux mode (its flux is
// overwritten by each interval's value and left at the last one); a
// constant-power system throws std::invalid_argument, since its matrix depends
// on the composition and no single A_k represents the interval.
//
// Every matrix carries a structural entry for every channel of the system --
// including channels with a zero cross section, and every entry of a zero-flux
// interval -- stored as an explicit zero. The values are those of
// assemble(); the structure is what rateSensitivities() integrates over and
// what reactionSensitivities() requires, so matrices built any other way can
// leave a channel's derivative unavailable.
std::vector<Eigen::SparseMatrix<double>> intervalMatrices(DepletionSystem& system,
                                                          const std::vector<double>& flux);

// --- First-order sensitivities ---------------------------------------------

// The integrand n*_i(t) n_j(t) is a product of two trajectories, so its
// integral is not a single matrix exponential; it is quadrature over the
// trajectories, with the forward and adjoint values at every node obtained by
// exact solves from the interval end points.
//
// Each schedule interval is split into `subIntervals` equal pieces and
// Gauss-Legendre with `gaussPoints` nodes (1..5) is applied to each. That
// alone converges only linearly, because both trajectories carry boundary
// layers at the interval ends: when the flux switches, every short-lived
// nuclide re-equilibrates over its own removal time (Xe-135 under flux: ~20
// minutes), and the adjoint does the same backward from the response weight.
// `endRefinements` resolves them by bisecting the first and last piece of each
// interval that many times toward the end, so the smallest piece is
// h / 2^endRefinements with h the uniform piece length. Raise it until the
// shortest removal time in the problem is no smaller than that; raise
// subIntervals for slow variation in the interior. The result should be
// checked against a finer rule once per problem class; the defaults suit a
// thermal pin with a few intervals of days to months.
struct SensitivityOptions {
  int gaussPoints = 5;
  int subIntervals = 4;
  int endRefinements = 6;
};

// S_k(i, j) = integral over interval k of n*_i(t) n_j(t) dt, one matrix per
// interval, restricted to the union of the A_k sparsity patterns -- every entry
// a decay constant or cross section can move. An entry outside that pattern is
// not an integral that came out zero, it is one that was never taken, and the
// contractions below refuse to report a derivative that depends on one; build
// A with intervalMatrices(), which carries the structure of every channel
// whatever its cross section. An empty schedule gives an empty result.
// dR/dn0 is not computed by this function: it is adj.nStar.front(). `fwd` and
// `adj` must come from depleteLinear() and
// adjointDeplete()/adjointDepleteIntegrated() over the same A and dts.
std::vector<Eigen::SparseMatrix<double>> rateSensitivities(
    const DepletionResult& fwd, const AdjointResult& adj,
    const std::vector<Eigen::SparseMatrix<double>>& A, const std::vector<double>& dts,
    SensitivityOptions options = {}, CramOrder order = CramOrder::CRAM48);

// dR/d(lambda_i) for every nuclide, indexed like the chain (0 for a nuclide
// with no decay). Contracts the summed S_k with dA/d(lambda_i): -1 on the
// diagonal, +branching at each mode's daughter, +branching*yield at each
// spontaneous-fission product.
Eigen::VectorXd decayConstantSensitivities(const DepletionChain& chain,
                                           const std::vector<Eigen::SparseMatrix<double>>& S);

// dR/d(sigma) [per barn] for every channel of the system, in the order of
// system.reactions(). Contracts each S_k with flux[k] * 1e-24 * dA/d(sigma):
// -1 on the parent diagonal, +1 at a tracked target, +yield at each fission
// product (yields at the channel's energy). flux.size() must equal S.size().
// Throws std::invalid_argument if S does not carry an entry for every one of
// those positions -- see rateSensitivities() above.
struct ReactionSensitivity {
  Zai parent;
  ReactionType type;
  std::optional<Zai> target;
  double dRdSigma;
};

std::vector<ReactionSensitivity> reactionSensitivities(
    const DepletionSystem& system, const std::vector<Eigen::SparseMatrix<double>>& S,
    const std::vector<double>& flux);

}  // namespace cram
