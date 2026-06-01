#pragma once
//
// Time integrators that march a depletion problem dn/dt = A(n) n over a burnup
// history. The burnup matrix A is composition-dependent in general (constant-
// power flux normalization makes the reaction rates depend on the current
// inventory), so a step cannot just reuse one exp(A*dt): the higher-order
// schemes re-evaluate A at predicted intermediate compositions.
//
// The algorithms and their coefficients mirror OpenMC's openmc.deplete
// integrators (Yu & Forget, Ann. Nucl. Energy 170 (2022) 108973, Table 1):
//   * Predictor          O(1)  -- one matrix exponential, A frozen at BOS
//   * CE/CM              O(2)  -- constant extrapolate to midpoint, then correct
//   * CE/LI (CFQ4)       O(2)  -- constant extrapolate, commutator-free linear
//   * LE/QI (CFQ4)       O(2)  -- linear extrapolate, commutator-free quadratic
//   * CF4                O(4)  -- commutator-free 4th-order Lie integrator
//
// CE/LI and LE/QI use the two-exponential CFQ4 forms used by OpenMC, not the
// single averaged-matrix forms; the matrix-combination coefficients and the
// order in which the two exponentials are applied match the OpenMC source so
// that, given identical A(n), this engine reproduces OpenMC's trajectory.
//
#include <Eigen/SparseCore>
#include <functional>
#include <memory>
#include <vector>

#include "cram/cram.hpp"

namespace cram {

// Assemble the burnup matrix A (per unit time, 1/s) for a given composition n.
// Reused at every intermediate composition a scheme evaluates.
using MatrixBuilder = std::function<Eigen::SparseMatrix<double>(const Eigen::VectorXd&)>;

enum class IntegratorKind { Predictor, CECM, CELI, LEQI, CF4 };

// One depletion integrator. step() advances n across a single interval of size
// dt [s] and may carry cross-step history (LE/QI needs the previous interval).
// Call reset() before starting a fresh march to clear that history.
class Integrator {
public:
  Integrator(MatrixBuilder build, CramOrder order) : build_(std::move(build)), order_(order) {}
  virtual ~Integrator() = default;

  virtual Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) = 0;
  virtual void reset() {}

  CramOrder order() const { return order_; }

protected:
  // exp(M * dt) * v, evaluated with CRAM. M is a per-unit-time matrix (or a
  // linear combination of such matrices); dt carries the interval length.
  Eigen::VectorXd expm(const Eigen::SparseMatrix<double>& M, double dt,
                       const Eigen::VectorXd& v) const {
    return cramSolve(M, v, dt, order_);
  }

  MatrixBuilder build_;
  CramOrder order_;
};

// Construct an integrator of the requested kind.
std::unique_ptr<Integrator> makeIntegrator(IntegratorKind kind, MatrixBuilder build,
                                           CramOrder order = CramOrder::CRAM48);

// Parse an integrator name ("predictor", "cecm", "celi", "leqi", "cf4";
// case-insensitive, '/' and '-' ignored). Throws std::invalid_argument if the
// name is not recognized.
IntegratorKind integratorKindFromName(const std::string& name);
const char* integratorName(IntegratorKind kind);

// --- Marching a step schedule ----------------------------------------------

// Per-point depletion trajectory: time[k] is the cumulative time (s) at point
// k and n[k] the composition there. Sizes are dts.size() + 1 (the initial
// point plus one per interval).
struct DepletionResult {
  std::vector<double> time;
  std::vector<Eigen::VectorXd> n;
};

// March n0 through the given intervals dts [s] using integ, recording the
// composition at every point. Calls integ.reset() first.
DepletionResult deplete(Integrator& integ, const Eigen::VectorXd& n0,
                        const std::vector<double>& dts);

}  // namespace cram
