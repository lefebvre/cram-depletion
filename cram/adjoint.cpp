#include "cram/adjoint.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

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

Eigen::SparseMatrix<double> sourceGenerator(const Eigen::SparseMatrix<double>& A,
                                            const Eigen::VectorXd& s) {
  if (A.rows() != A.cols() || s.size() != A.rows())
    throw std::invalid_argument("cram: sourceGenerator needs a square A and a source of size N");
  const Eigen::Index N = A.rows();
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(A.nonZeros() + N));
  for (Eigen::Index k = 0; k < A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
      t.emplace_back(it.row(), it.col(), it.value());
  for (Eigen::Index i = 0; i < N; ++i)
    if (s(i) != 0.0)
      t.emplace_back(i, N, s(i));
  Eigen::SparseMatrix<double> G(N + 1, N + 1);
  G.setFromTriplets(t.begin(), t.end());
  G.makeCompressed();
  return G;
}

Eigen::SparseMatrix<double> augmentedGenerator(const Eigen::SparseMatrix<double>& A) {
  if (A.rows() != A.cols())
    throw std::invalid_argument("cram: augmentedGenerator needs a square A");
  const Eigen::Index N = A.rows();
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(A.nonZeros() + N));
  for (Eigen::Index k = 0; k < A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
      t.emplace_back(it.row(), it.col(), it.value());
  for (Eigen::Index i = 0; i < N; ++i)
    t.emplace_back(N + i, i, 1.0);
  Eigen::SparseMatrix<double> G(2 * N, 2 * N);
  G.setFromTriplets(t.begin(), t.end());
  G.makeCompressed();
  return G;
}

namespace {

void requireSchedule(const std::vector<Eigen::SparseMatrix<double>>& A,
                     const std::vector<double>& dts) {
  if (A.size() != dts.size())
    throw std::invalid_argument("cram: one matrix per interval is required (" +
                                std::to_string(A.size()) + " matrices, " +
                                std::to_string(dts.size()) + " intervals)");
}

// One step of an inhomogeneous problem: exp(G dt) applied to (v; 1) for
// G = sourceGenerator(M, s), with the trailing 1 dropped, which is
// exp(M dt) v + integral_0^dt exp(M u) s du.
Eigen::VectorXd stepWithSource(CramSolver& solver, const Eigen::SparseMatrix<double>& G,
                               const Eigen::VectorXd& v, double dt) {
  const Eigen::Index N = G.rows() - 1;
  solver.prepare(G, dt);
  Eigen::VectorXd x(N + 1);
  x.head(N) = v;
  x(N) = 1.0;
  return solver.apply(x).head(N);
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

Eigen::VectorXd integratedInventory(const std::vector<Eigen::SparseMatrix<double>>& A,
                                    const std::vector<double>& dts, const Eigen::VectorXd& n0,
                                    CramOrder order) {
  requireSchedule(A, dts);
  const Eigen::Index N = n0.size();
  CramSolver solver(order);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(2 * N);
  x.head(N) = n0;
  for (std::size_t k = 0; k < dts.size(); ++k) {
    // The bottom block accumulates across intervals: exp(G dt) leaves it
    // untouched apart from adding this interval's integral.
    solver.prepare(augmentedGenerator(A[k]), dts[k]);
    x = solver.apply(x);
  }
  return x.tail(N);
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

AdjointResult adjointDepleteIntegrated(const std::vector<Eigen::SparseMatrix<double>>& A,
                                       const std::vector<double>& dts, const Eigen::VectorXd& w,
                                       CramOrder order) {
  requireSchedule(A, dts);
  const std::size_t K = dts.size();

  AdjointResult r;
  r.time.assign(K + 1, 0.0);
  for (std::size_t k = 0; k < K; ++k)
    r.time[k + 1] = r.time[k] + dts[k];
  r.nStar.assign(K + 1, Eigen::VectorXd());
  r.nStar[K] = Eigen::VectorXd::Zero(w.size());
  r.integratedWeight = w;

  // g(t_k) = exp(A_k^T dt_k) g(t_{k+1}) + integral_0^{dt_k} exp(A_k^T u) w du,
  // which is one solve on the source generator of (A_k^T, w).
  CramSolver solver(order);
  for (std::size_t k = K; k-- > 0;)
    r.nStar[k] =
        stepWithSource(solver, sourceGenerator(transposed(A[k]), w), r.nStar[k + 1], dts[k]);
  return r;
}

namespace {

// Every entry of A that some channel of `system` can move, as explicit zeros:
// the parent diagonal, a tracked target, and each product of the yield table a
// fission channel would use.
//
// assemble() emits nothing at all for a channel whose rate is zero -- a cross
// section of zero, or a zero-flux interval -- so those entries would be absent
// from the sparsity pattern, and rateSensitivities() integrates only the
// pattern. An absent entry and a zero one are the same to coeff() but mean
// opposite things: an integral that came out zero, versus one that was never
// taken. Carrying the structure keeps the two apart, at the price of a handful
// of stored zeros in a matrix whose values are unchanged.
std::vector<Eigen::Triplet<double>> channelPattern(const DepletionSystem& system) {
  const DepletionChain& chain = system.chain();
  std::vector<Eigen::Triplet<double>> t;
  std::vector<std::pair<int, double>> column;
  for (const auto& [parent, rxns] : system.reactions()) {
    const int p = chain.indexOf(parent);
    for (const ReactionXS& r : rxns) {
      // The entries assemble() writes for this channel whenever its rate is
      // not zero -- the same stencil, so the two cannot disagree.
      chain.reactionStencil(parent, r.type == ReactionType::Fission, r.target, r.energy, column);
      for (const auto& [j, c] : column)
        t.emplace_back(j, p, 0.0);
    }
  }
  return t;
}

// A with `extra` merged into its structure. setFromTriplets() sums duplicates
// and stores what it is given, zero or not, so the values are exactly A's.
Eigen::SparseMatrix<double> withPattern(const Eigen::SparseMatrix<double>& A,
                                        const std::vector<Eigen::Triplet<double>>& extra) {
  std::vector<Eigen::Triplet<double>> t = extra;
  t.reserve(extra.size() + static_cast<std::size_t>(A.nonZeros()));
  for (Eigen::Index k = 0; k < A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
      t.emplace_back(it.row(), it.col(), it.value());
  Eigen::SparseMatrix<double> out(A.rows(), A.cols());
  out.setFromTriplets(t.begin(), t.end());
  out.makeCompressed();
  return out;
}

}  // namespace

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
  const std::vector<Eigen::Triplet<double>> channels = channelPattern(system);
  // Each interval is assembled at its own flux, and the caller's flux is put
  // back afterwards -- on the error path too -- so building the schedule does
  // not leave the system describing its last interval.
  const double callerFlux = system.fluxFor(any);
  try {
    for (double phi : flux) {
      system.setConstantFlux(phi);
      A.push_back(withPattern(system.assemble(any), channels));
    }
  } catch (...) {
    system.setConstantFlux(callerFlux);
    throw;
  }
  system.setConstantFlux(callerFlux);
  return A;
}

// --- Sensitivities ------------------------------------------------------------

namespace {

struct GaussRule {
  std::vector<double> x;  // nodes on [-1, 1]
  std::vector<double> w;  // weights summing to 2
};

GaussRule gaussLegendre(int points) {
  switch (points) {
    case 1:
      return {{0.0}, {2.0}};
    case 2:
      return {{-0.5773502691896257, 0.5773502691896257}, {1.0, 1.0}};
    case 3:
      return {{-0.7745966692414834, 0.0, 0.7745966692414834},
              {0.5555555555555556, 0.8888888888888888, 0.5555555555555556}};
    case 4:
      return {{-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526},
              {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538}};
    case 5:
      return {
          {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640},
          {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
           0.2369268850561891}};
    default:
      throw std::invalid_argument("cram: SensitivityOptions::gaussPoints must be 1..5");
  }
}

// The quadrature pieces of one interval [0, dt]: `sub` uniform pieces, with
// the first and last bisected `levels` times toward the interval end so a
// boundary layer there is resolved geometrically.
std::vector<std::pair<double, double>> quadraturePieces(double dt, int sub, int levels) {
  std::vector<std::pair<double, double>> pieces;
  // A single piece is graded at both ends by treating it as two halves.
  const int n = (sub == 1 && levels > 0) ? 2 : sub;
  const double h = dt / n;
  pieces.reserve(static_cast<std::size_t>(n));
  for (int p = 0; p < n; ++p)
    pieces.emplace_back(p * h, (p + 1) * h);
  if (levels == 0)
    return pieces;

  std::vector<std::pair<double, double>> graded;
  // First piece [0, h]: [0, h/2^L], [h/2^L, h/2^(L-1)], ..., [h/2, h].
  double lo = 0.0;
  for (int m = levels; m >= 1; --m) {
    const double hi = h / static_cast<double>(1 << m);
    graded.emplace_back(lo, hi);
    lo = hi;
  }
  graded.emplace_back(lo, h);
  for (int p = 1; p + 1 < n; ++p)
    graded.push_back(pieces[static_cast<std::size_t>(p)]);
  // Last piece [dt-h, dt], mirrored: [dt-h, dt-h/2], ..., [dt-h/2^L, dt].
  double hi = dt;
  std::vector<std::pair<double, double>> tail;
  for (int m = levels; m >= 1; --m) {
    const double lo2 = dt - h / static_cast<double>(1 << m);
    tail.emplace_back(lo2, hi);
    hi = lo2;
  }
  graded.emplace_back(dt - h, hi);
  for (auto it = tail.rbegin(); it != tail.rend(); ++it)
    graded.push_back(*it);
  return graded;
}

// The union of the sparsity patterns, as a matrix of zeros with that structure.
// Empty in, empty out: an empty schedule has no dimension to report.
Eigen::SparseMatrix<double> patternUnion(const std::vector<Eigen::SparseMatrix<double>>& A) {
  if (A.empty())
    return {};
  std::vector<Eigen::Triplet<double>> t;
  for (const auto& M : A)
    for (Eigen::Index k = 0; k < M.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
        t.emplace_back(it.row(), it.col(), 0.0);
  Eigen::SparseMatrix<double> P(A.front().rows(), A.front().cols());
  P.setFromTriplets(t.begin(), t.end());
  P.makeCompressed();
  return P;
}

}  // namespace

std::vector<Eigen::SparseMatrix<double>> rateSensitivities(
    const DepletionResult& fwd, const AdjointResult& adj,
    const std::vector<Eigen::SparseMatrix<double>>& A, const std::vector<double>& dts,
    SensitivityOptions options, CramOrder order) {
  requireSchedule(A, dts);
  if (fwd.n.size() != dts.size() + 1 || adj.nStar.size() != dts.size() + 1)
    throw std::invalid_argument(
        "cram: rateSensitivities needs trajectories with one point per schedule point");
  if (options.subIntervals < 1)
    throw std::invalid_argument("cram: SensitivityOptions::subIntervals must be >= 1");
  if (options.endRefinements < 0 || options.endRefinements > 30)
    throw std::invalid_argument("cram: SensitivityOptions::endRefinements must be 0..30");
  const GaussRule rule = gaussLegendre(options.gaussPoints);
  // An empty schedule is a valid one everywhere else here -- depleteLinear()
  // and adjointDeplete() return their single initial point for it -- so it
  // returns no interval matrices rather than reaching for A.front() below.
  if (dts.empty())
    return {};
  const Eigen::SparseMatrix<double> pattern = patternUnion(A);

  std::vector<Eigen::SparseMatrix<double>> S;
  S.reserve(dts.size());
  const Eigen::Index N = fwd.n.front().size();
  const bool integrated = adj.integratedWeight.has_value();
  CramSolver solver(order);
  std::vector<std::pair<double, double>> nodes;  // (offset from interval start, weight)
  for (std::size_t k = 0; k < dts.size(); ++k) {
    Eigen::SparseMatrix<double> Sk = pattern;  // zeros with the union structure

    // Each node needs n at offset tau, exp(A tau) n_k, and n* at offset tau,
    // exp(A^T (dt - tau)) n*_{k+1}. The node set is symmetric about the
    // interval midpoint -- the pieces mirror each other and so do the Gauss
    // points -- so dt - tau is the offset of the mirror node, and the transpose
    // of the factorization prepared for exp(A tau_mirror) gives the adjoint at
    // this one. One prepare() per node therefore serves a forward value and an
    // adjoint value, half of what a forward and a backward solver each paying
    // for their own would cost, with no marching error.
    //
    // For a time-integrated response the adjoint carries a source:
    // n*(tau) = [exp(G (dt - tau)) (n*_{k+1}; 1)].head(N) with
    // G = sourceGenerator(A^T, w) = [[A^T, w], [0, 0]]. Its transpose
    // [[A, 0], [w^T, 0]] applied to (n_k; 0) yields exp(A tau) n_k in the top
    // block, so the same pairing works on that (N+1)-sized matrix.
    const Eigen::SparseMatrix<double> M =
        integrated ? transposed(sourceGenerator(transposed(A[k]), *adj.integratedWeight)) : A[k];
    Eigen::VectorXd forwardStart = Eigen::VectorXd::Zero(M.rows());
    Eigen::VectorXd adjointEnd = Eigen::VectorXd::Zero(M.rows());
    forwardStart.head(N) = fwd.n[k];
    adjointEnd.head(N) = adj.nStar[k + 1];
    if (integrated)
      adjointEnd(N) = 1.0;

    nodes.clear();
    for (const auto& [a, b] :
         quadraturePieces(dts[k], options.subIntervals, options.endRefinements)) {
      const double h = b - a;
      for (std::size_t g = 0; g < rule.x.size(); ++g)
        nodes.emplace_back(a + 0.5 * (1.0 + rule.x[g]) * h, 0.5 * h * rule.w[g]);
    }

    const auto accumulate = [&Sk](double weight, const Eigen::VectorXd& nStar,
                                  const Eigen::VectorXd& n) {
      for (Eigen::Index c = 0; c < Sk.outerSize(); ++c)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Sk, c); it; ++it)
          it.valueRef() += weight * nStar(it.row()) * n(it.col());
    };
    const std::size_t m = nodes.size();
    for (std::size_t g = 0; g < (m + 1) / 2; ++g) {
      const std::size_t mirror = m - 1 - g;
      const auto [tau, weight] = nodes[g];
      const auto [tauMirror, weightMirror] = nodes[mirror];
      if (std::abs(tauMirror - (dts[k] - tau)) > 1e-12 * dts[k])
        throw std::logic_error("cram: rateSensitivities quadrature nodes are not symmetric");

      solver.prepare(M, tau);
      const Eigen::VectorXd n = solver.apply(forwardStart).head(N);
      const Eigen::VectorXd nStarMirror = solver.applyTranspose(adjointEnd).head(N);
      if (mirror == g) {  // the midpoint node is its own mirror
        accumulate(weight, nStarMirror, n);
        continue;
      }
      solver.prepare(M, tauMirror);
      const Eigen::VectorXd nMirror = solver.apply(forwardStart).head(N);
      const Eigen::VectorXd nStar = solver.applyTranspose(adjointEnd).head(N);
      accumulate(weight, nStar, n);
      accumulate(weightMirror, nStarMirror, nMirror);
    }
    S.push_back(std::move(Sk));
  }
  return S;
}

namespace {

// A stencil with its repeated rows merged, so that a product which is also the
// parent nets against the removal instead of replacing it, and so the
// contraction below can be a single walk of one column.
std::unordered_map<int, double> mergeStencil(const std::vector<std::pair<int, double>>& stencil) {
  std::unordered_map<int, double> gain;
  for (const auto& [row, c] : stencil)
    gain[row] += c;
  return gain;
}

// sum over `gain` of coefficient * S(row, col), walking column `col` of S
// once. Empty when S carries no entry for one of the rows. Reading them with
// coeff() instead would report 0 for an entry S never integrated --
// indistinguishable from one that integrated to zero -- and, since the removal
// term is present either way, would hand back part of a derivative with no
// sign left to trust.
std::optional<double> contractColumn(const Eigen::SparseMatrix<double>& S, int col,
                                     const std::unordered_map<int, double>& gain) {
  double dA = 0.0;
  std::size_t found = 0;
  for (Eigen::SparseMatrix<double>::InnerIterator it(S, col); it; ++it) {
    const auto g = gain.find(static_cast<int>(it.row()));
    if (g == gain.end())
      continue;
    dA += g->second * it.value();
    ++found;
  }
  if (found != gain.size())
    return std::nullopt;
  return dA;
}

void requireChainSized(const std::vector<Eigen::SparseMatrix<double>>& S,
                       const DepletionChain& chain, const char* who) {
  const bool wrongSize =
      std::any_of(S.begin(), S.end(), [&chain](const Eigen::SparseMatrix<double>& Sk) {
        return Sk.rows() != chain.size() || Sk.cols() != chain.size();
      });
  if (wrongSize)
    throw std::invalid_argument(std::string("cram: ") + who +
                                " needs each S_k sized like the chain");
}

constexpr const char* kBuildWithIntervalMatrices =
    "S_k carries only the sparsity of the matrices it was built from; build them with "
    "intervalMatrices(), which keeps a structural entry for every decay and every channel "
    "of the system";

}  // namespace

Eigen::VectorXd decayConstantSensitivities(const DepletionChain& chain,
                                           const std::vector<Eigen::SparseMatrix<double>>& S) {
  requireChainSized(S, chain, "decayConstantSensitivities");
  Eigen::VectorXd out = Eigen::VectorXd::Zero(chain.size());
  if (S.empty())
    return out;
  Eigen::SparseMatrix<double> total = S.front();
  for (std::size_t k = 1; k < S.size(); ++k)
    total = total + S[k];

  // dA/d(lambda_i) is exactly the column decayMatrix() scales by lambda_i.
  std::vector<std::pair<int, double>> column;
  for (int i = 0; i < chain.size(); ++i) {
    chain.decayStencil(i, column);
    if (column.empty())
      continue;  // no decay: lambda_i appears nowhere in A
    const std::optional<double> s = contractColumn(total, i, mergeStencil(column));
    if (!s)
      throw std::invalid_argument(
          "cram: decayConstantSensitivities: S has no entry for part of the decay of " +
          chain.nuclides()[static_cast<std::size_t>(i)].str() +
          ", so its derivative would be incomplete. " + kBuildWithIntervalMatrices);
    out(i) = *s;
  }
  return out;
}

std::vector<ReactionSensitivity> reactionSensitivities(
    const DepletionSystem& system, const std::vector<Eigen::SparseMatrix<double>>& S,
    const std::vector<double>& flux) {
  if (S.size() != flux.size())
    throw std::invalid_argument("cram: reactionSensitivities needs one flux per S_k");
  constexpr double kBarn = 1e-24;
  const DepletionChain& chain = system.chain();
  requireChainSized(S, chain, "reactionSensitivities");

  std::vector<ReactionSensitivity> out;
  std::vector<std::pair<int, double>> column;
  for (const auto& [parent, rxns] : system.reactions()) {
    const int p = chain.indexOf(parent);
    for (const ReactionXS& r : rxns) {
      double s = 0.0;
      if (p >= 0) {
        // dA/dsigma is flux * 1e-24 times the stencil assemble() scales by the
        // reaction rate.
        chain.reactionStencil(parent, r.type == ReactionType::Fission, r.target, r.energy, column);
        const std::unordered_map<int, double> gain = mergeStencil(column);
        for (std::size_t k = 0; k < S.size(); ++k) {
          const std::optional<double> dA = contractColumn(S[k], p, gain);
          if (!dA)
            throw std::invalid_argument(
                "cram: reactionSensitivities: S has no entry for part of the " +
                std::string(reactionName(r.type)) + " channel of " + parent.str() +
                ", so its derivative would be incomplete; a channel with no rate on any "
                "interval contributes no structure. " +
                kBuildWithIntervalMatrices);
          s += flux[k] * kBarn * *dA;
        }
      }
      out.push_back(
          ReactionSensitivity{.parent = parent, .type = r.type, .target = r.target, .dRdSigma = s});
    }
  }
  return out;
}

}  // namespace cram
