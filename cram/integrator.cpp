#include "cram/integrator.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace cram {
namespace {

using SpMat = Eigen::SparseMatrix<double>;

// First-order predictor (OpenMC PredictorIntegrator):
//   n_{i+1} = exp(dt * A(n_i)) n_i
class PredictorIntegrator final : public Integrator {
public:
  using Integrator::Integrator;
  Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) override {
    return expm(build_(n), dt, n);
  }
};

// CE/CM predictor-corrector (OpenMC CECMIntegrator):
//   n_{1/2} = exp(dt/2 * A(n_i)) n_i
//   n_{i+1} = exp(dt   * A(n_{1/2})) n_i
class CECMIntegrator final : public Integrator {
public:
  using Integrator::Integrator;
  Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) override {
    SpMat A0 = build_(n);
    Eigen::VectorXd nMid = expm(A0, dt / 2.0, n);
    return expm(build_(nMid), dt, n);
  }
};

// CE/LI CFQ4 predictor-corrector (OpenMC CELIIntegrator). Predictor is constant
// extrapolation; corrector is the two-exponential commutator-free linear form.
// Application order follows the OpenMC source: f1-combination first, then f2.
//   n^p     = exp(dt * A(n_i)) n_i
//   n_inter = exp(dt * (5/12 A(n_i) + 1/12 A(n^p))) n_i
//   n_{i+1} = exp(dt * (1/12 A(n_i) + 5/12 A(n^p))) n_inter
class CELIIntegrator final : public Integrator {
public:
  using Integrator::Integrator;
  Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) override {
    SpMat A0 = build_(n);
    Eigen::VectorXd nCe = expm(A0, dt, n);  // predictor
    SpMat A1 = build_(nCe);
    SpMat f1 = (5.0 / 12.0) * A0 + (1.0 / 12.0) * A1;
    Eigen::VectorXd nInter = expm(f1, dt, n);
    SpMat f2 = (1.0 / 12.0) * A0 + (5.0 / 12.0) * A1;
    return expm(f2, dt, nInter);
  }
};

// CF4 commutator-free 4th-order Lie integrator (OpenMC CF4Integrator). Five
// matrix exponentials; A_k below denotes dt*A evaluated at the kth intermediate
// composition. Application order follows the OpenMC source.
class CF4Integrator final : public Integrator {
public:
  using Integrator::Integrator;
  Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) override {
    SpMat A0 = build_(n);
    Eigen::VectorXd n1 = expm(0.5 * A0, dt, n);  // exp(A1/2) n
    SpMat A1 = build_(n1);
    Eigen::VectorXd n2 = expm(0.5 * A1, dt, n);  // exp(A2/2) n
    SpMat A2 = build_(n2);
    SpMat f2 = -0.5 * A0 + A2;              // -A1/2 + A3
    Eigen::VectorXd n3 = expm(f2, dt, n1);  // exp(-A1/2 + A3) n1
    SpMat A3 = build_(n3);
    SpMat f3 = 0.25 * A0 + (1.0 / 6.0) * A1 + (1.0 / 6.0) * A2 - (1.0 / 12.0) * A3;
    Eigen::VectorXd nInter = expm(f3, dt, n);
    SpMat f4 = -(1.0 / 12.0) * A0 + (1.0 / 6.0) * A1 + (1.0 / 6.0) * A2 + 0.25 * A3;
    return expm(f4, dt, nInter);
  }
};

// LE/QI CFQ4 predictor-corrector (OpenMC LEQIIntegrator). Needs the burnup
// matrix from the previous interval's start-of-step composition (A_{-1}) and
// the previous interval length (dt_l). The first step (no history) falls back
// to CE/LI, matching OpenMC. Application order follows the OpenMC source.
class LEQIIntegrator final : public Integrator {
public:
  using Integrator::Integrator;

  void reset() override { hasPrev_ = false; }

  Eigen::VectorXd step(const Eigen::VectorXd& n, double dt) override {
    SpMat A0 = build_(n);  // A(n_i)

    if (!hasPrev_) {
      // CE/LI bootstrap, but remember this step's BOS matrix for the next one.
      Eigen::VectorXd nCe = expm(A0, dt, n);
      SpMat A1 = build_(nCe);
      Eigen::VectorXd nInter = expm((5.0 / 12.0) * A0 + (1.0 / 12.0) * A1, dt, n);
      Eigen::VectorXd nEnd = expm((1.0 / 12.0) * A0 + (5.0 / 12.0) * A1, dt, nInter);
      prevA0_ = A0;
      prevDt_ = dt;
      hasPrev_ = true;
      return nEnd;
    }

    const double dtl = prevDt_;
    const SpMat& Am1 = prevA0_;  // A(n_{i-1})

    // Predictor (linear extrapolation, two-exponential CFQ4 form).
    SpMat f1 = (-dt / (12.0 * dtl)) * Am1 + ((dt + 6.0 * dtl) / (12.0 * dtl)) * A0;
    Eigen::VectorXd nInter = expm(f1, dt, n);
    SpMat f2 = (-5.0 * dt / (12.0 * dtl)) * Am1 + ((5.0 * dt + 6.0 * dtl) / (12.0 * dtl)) * A0;
    Eigen::VectorXd nPred = expm(f2, dt, nInter);

    SpMat A1 = build_(nPred);  // A(n^p_{i+1})

    // Corrector (quadratic interpolation, two-exponential CFQ4 form).
    const double den = 12.0 * dtl * (dtl + dt);
    SpMat f3 = (-dt * dt / den) * Am1 + ((dt * dt + 6.0 * dt * dtl + 5.0 * dtl * dtl) / den) * A0 +
               (dtl / (12.0 * (dtl + dt))) * A1;
    Eigen::VectorXd nInter2 = expm(f3, dt, n);
    SpMat f4 = (-dt * dt / den) * Am1 + ((dt * dt + 2.0 * dt * dtl + dtl * dtl) / den) * A0 +
               ((4.0 * dt + 5.0 * dtl) / (12.0 * (dtl + dt))) * A1;
    Eigen::VectorXd nEnd = expm(f4, dt, nInter2);

    prevA0_ = A0;
    prevDt_ = dt;
    return nEnd;
  }

private:
  bool hasPrev_ = false;
  double prevDt_ = 0.0;
  SpMat prevA0_;
};

}  // namespace

std::unique_ptr<Integrator> makeIntegrator(IntegratorKind kind, MatrixBuilder build,
                                           CramOrder order) {
  switch (kind) {
    case IntegratorKind::Predictor:
      return std::make_unique<PredictorIntegrator>(std::move(build), order);
    case IntegratorKind::CECM:
      return std::make_unique<CECMIntegrator>(std::move(build), order);
    case IntegratorKind::CELI:
      return std::make_unique<CELIIntegrator>(std::move(build), order);
    case IntegratorKind::LEQI:
      return std::make_unique<LEQIIntegrator>(std::move(build), order);
    case IntegratorKind::CF4:
      return std::make_unique<CF4Integrator>(std::move(build), order);
  }
  throw std::invalid_argument("cram: unknown IntegratorKind");
}

IntegratorKind integratorKindFromName(const std::string& name) {
  std::string key;
  for (char c : name) {
    if (c == '/' || c == '-' || c == '_' || c == ' ')
      continue;
    key.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  }
  if (key == "predictor" || key == "ce")
    return IntegratorKind::Predictor;
  if (key == "cecm")
    return IntegratorKind::CECM;
  if (key == "celi")
    return IntegratorKind::CELI;
  if (key == "leqi")
    return IntegratorKind::LEQI;
  if (key == "cf4")
    return IntegratorKind::CF4;
  throw std::invalid_argument("cram: unknown integrator name '" + name + "'");
}

const char* integratorName(IntegratorKind kind) {
  switch (kind) {
    case IntegratorKind::Predictor:
      return "predictor";
    case IntegratorKind::CECM:
      return "cecm";
    case IntegratorKind::CELI:
      return "celi";
    case IntegratorKind::LEQI:
      return "leqi";
    case IntegratorKind::CF4:
      return "cf4";
  }
  return "unknown";
}

DepletionResult deplete(Integrator& integ, const Eigen::VectorXd& n0,
                        const std::vector<double>& dts) {
  integ.reset();
  DepletionResult r;
  r.time.reserve(dts.size() + 1);
  r.n.reserve(dts.size() + 1);

  double t = 0.0;
  Eigen::VectorXd n = n0;
  r.time.push_back(t);
  r.n.push_back(n);
  for (double dt : dts) {
    n = integ.step(n, dt);
    t += dt;
    r.time.push_back(t);
    r.n.push_back(n);
  }
  return r;
}

}  // namespace cram
