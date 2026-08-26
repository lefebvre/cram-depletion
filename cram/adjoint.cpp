#include "cram/adjoint.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>

#include "cram/cram_solver.hpp"

namespace cram {

Eigen::SparseMatrix<double> transposed(const Eigen::SparseMatrix<double>& A) {
  Eigen::SparseMatrix<double> At = A.transpose();
  At.makeCompressed();
  return At;
}

Eigen::VectorXd cramSolveAdjoint(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& w,
                                 double dt, CramOrder order) {
  return cramSolve(transposed(A), w, dt, order);
}

namespace {

void requireSchedule(const std::vector<Eigen::SparseMatrix<double>>& A,
                     const std::vector<double>& dts) {
  if (A.size() != dts.size())
    throw std::invalid_argument("cram: one matrix per interval is required (" +
                                std::to_string(A.size()) + " matrices, " +
                                std::to_string(dts.size()) + " intervals)");
}

}  // namespace

DepletionResult depleteLinear(const std::vector<Eigen::SparseMatrix<double>>& A,
                              const std::vector<double>& dts, const Eigen::VectorXd& n0,
                              CramOrder order) {
  requireSchedule(A, dts);
  DepletionResult r;
  r.time.reserve(dts.size() + 1);
  r.n.reserve(dts.size() + 1);
  r.time.push_back(0.0);
  r.n.push_back(n0);

  CramSolver solver(order);
  double t = 0.0;
  for (std::size_t k = 0; k < dts.size(); ++k) {
    solver.prepare(A[k], dts[k]);
    r.n.push_back(solver.apply(r.n.back()));
    t += dts[k];
    r.time.push_back(t);
  }
  return r;
}

AdjointResult adjointDeplete(const std::vector<Eigen::SparseMatrix<double>>& A,
                             const std::vector<double>& dts, const Eigen::VectorXd& w,
                             CramOrder order) {
  requireSchedule(A, dts);
  const std::size_t K = dts.size();

  AdjointResult r;
  r.time.assign(K + 1, 0.0);
  for (std::size_t k = 0; k < K; ++k)
    r.time[k + 1] = r.time[k] + dts[k];
  r.nStar.assign(K + 1, Eigen::VectorXd());
  r.nStar[K] = w;

  // Backward in time: the transposed matrices share one sparsity pattern, so
  // every interval after the first reuses the symbolic analysis.
  CramSolver solver(order);
  for (std::size_t k = K; k-- > 0;) {
    solver.prepare(transposed(A[k]), dts[k]);
    r.nStar[k] = solver.apply(r.nStar[k + 1]);
  }
  return r;
}

std::vector<Eigen::SparseMatrix<double>> intervalMatrices(DepletionSystem& system,
                                                          const std::vector<double>& flux) {
  if (system.normalization() != DepletionSystem::Normalization::ConstantFlux)
    throw std::invalid_argument(
        "cram: intervalMatrices() needs a constant-flux system; under constant power the "
        "matrix depends on the composition and no single matrix represents an interval");
  std::vector<Eigen::SparseMatrix<double>> A;
  A.reserve(flux.size());
  // The composition is irrelevant under constant flux; any vector of the right
  // size will do.
  const Eigen::VectorXd any = Eigen::VectorXd::Zero(system.chain().size());
  for (double phi : flux) {
    system.setConstantFlux(phi);
    A.push_back(system.assemble(any));
  }
  return A;
}

}  // namespace cram
