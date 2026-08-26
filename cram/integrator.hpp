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
#include <string>
#include <vector>

#include "cram/cram.hpp"
#include "cram/cram_solver.hpp"

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
  Integrator(MatrixBuilder build, CramOrder order)
      : build_(std::move(build)), order_(order), solver_(order) {}
  virtual ~Integrator() = default;

  virtual Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) = 0;
  virtual void reset() {}

  CramOrder order() const { return order_; }

protected:
  // exp(M * dt) * v, evaluated with CRAM. M is a per-unit-time matrix (or a
  // linear combination of such matrices); dt carries the interval length.
  //
  // Routed through a CramSolver this integrator owns rather than the one-shot
  // cramSolve(), because every matrix a march produces shares one sparsity
  // pattern -- the burnup matrix's, at whatever composition it was built, and
  // any linear combination of those -- so the symbolic analysis inside the pole
  // factorizations is needed once per march, not once per exponential. It is
  // ~40% of a prepare() at N~1700, and a CF4 step spends five exponentials, so
  // a 100-step march was re-deriving it 500 times. Measured on the burnup shape
  // at N=1675 with CRAM48, over a sequence of same-pattern matrices: 20.9 ms
  // per exponential through cramSolve(), 15.2 ms through the cached solver --
  // 27% off every exponential after the first. That first prepare() costs more
  // than one cramSolve() (it pays the analysis and caches K factorizations), so
  // the trade breaks even at ~10 exponentials -- two CF4 steps.
  //
  // What it costs is memory: a prepared solver holds K complex LU
  // factorizations (24 for CRAM48), a few MB at N~1700, for as long as the
  // integrator lives. A many-region sweep should therefore build its integrator
  // inside the worker that marches the region -- which it must anyway, since
  // the builder captures that region's system -- rather than holding one per
  // region alive at once.
  Eigen::VectorXd expm(const Eigen::SparseMatrix<double>& M, double dt, const Eigen::VectorXd& v) {
    // As cramSolve(): exp(A*0) = I and exp(0*dt) = I, and neither is worth a
    // factorization.
    if (dt == 0.0 || M.nonZeros() == 0)
      return v;
    solver_.prepare(M, dt);
    return solver_.apply(v);
  }

  // One CE/LI CFQ4 step, from a start-of-step matrix that has already been
  // built:
  //   n^p     = exp(dt * A0) n
  //   n_inter = exp(dt * (5/12 A0 + 1/12 A(n^p))) n
  //   n_end   = exp(dt * (1/12 A0 + 5/12 A(n^p))) n_inter
  //
  // Shared because LE/QI's first step IS a CE/LI step -- OpenMC bootstraps it
  // the same way -- and LE/QI has already built A0 for its own history. Two
  // copies of the CFQ4 coefficients would be two places to correct against the
  // OpenMC source, and the convergence test only observes the aggregate order,
  // so a drift between them would not show up as a failure.
  Eigen::VectorXd celiStep(const Eigen::SparseMatrix<double>& A0, const Eigen::VectorXd& n,
                           double dt);

  MatrixBuilder build_;
  CramOrder order_;

private:
  CramSolver solver_;
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
