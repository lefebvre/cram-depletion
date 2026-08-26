#include "cram/adjoint.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
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

// exp(G dt) applied to (v; 1) for G = sourceGenerator(M, s): the first N
// entries of the result are exp(M dt) v + integral_0^dt exp(M u) s du.
// One step of an inhomogeneous problem through an ALREADY BUILT (N+1)x(N+1)
// generator G = sourceGenerator(M, s): exp(G dt) applied to (v; 1), with the
// trailing 1 dropped.
//
// G is a parameter rather than built here because only `dt` varies across the
// calls that matter. rateSensitivities() steps from one interval's generator to
// every quadrature node in it -- 80 of them at the default rule, 200 at the
// settings the tests use -- and building it per node copied the whole triplet
// list of A each time for a matrix that never changed.
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
  for (const auto& [parent, rxns] : system.reactions()) {
    const int p = chain.indexOf(parent);
    if (p < 0)
      continue;
    for (const ReactionXS& r : rxns) {
      t.emplace_back(p, p, 0.0);  // removal of the parent
      if (r.type == ReactionType::Fission) {
        if (const FissionYields* y = chain.nearestYields(parent, r.energy)) {
          for (const auto& [product, yield] : y->products) {
            const int j = chain.indexOf(product);
            if (j >= 0)
              t.emplace_back(j, p, 0.0);
          }
        }
      } else if (r.target) {
        const int j = chain.indexOf(*r.target);
        if (j >= 0)
          t.emplace_back(j, p, 0.0);
      }
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
  for (double phi : flux) {
    system.setConstantFlux(phi);
    A.push_back(withPattern(system.assemble(any), channels));
  }
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
  CramSolver forward(order);
  CramSolver backward(order);
  for (std::size_t k = 0; k < dts.size(); ++k) {
    Eigen::SparseMatrix<double> Sk = pattern;  // zeros with the union structure
    const Eigen::SparseMatrix<double> At = transposed(A[k]);
    // Depends only on the interval, so it is built here and not at every node.
    const Eigen::SparseMatrix<double> source = adj.integratedWeight
                                                   ? sourceGenerator(At, *adj.integratedWeight)
                                                   : Eigen::SparseMatrix<double>();
    for (const auto& [a, b] :
         quadraturePieces(dts[k], options.subIntervals, options.endRefinements)) {
      const double h = b - a;
      for (std::size_t g = 0; g < rule.x.size(); ++g) {
        // Offset of this node from the interval start, and from its end.
        const double tau = a + 0.5 * (1.0 + rule.x[g]) * h;
        const double rest = dts[k] - tau;
        forward.prepare(A[k], tau);
        const Eigen::VectorXd n = forward.apply(fwd.n[k]);
        Eigen::VectorXd nStar;
        if (adj.integratedWeight) {
          nStar = stepWithSource(backward, source, adj.nStar[k + 1], rest);
        } else {
          backward.prepare(At, rest);
          nStar = backward.apply(adj.nStar[k + 1]);
        }
        const double weight = 0.5 * h * rule.w[g];
        for (Eigen::Index c = 0; c < Sk.outerSize(); ++c)
          for (Eigen::SparseMatrix<double>::InnerIterator it(Sk, c); it; ++it)
            it.valueRef() += weight * nStar(it.row()) * n(it.col());
      }
    }
    S.push_back(std::move(Sk));
  }
  return S;
}

Eigen::VectorXd decayConstantSensitivities(const DepletionChain& chain,
                                           const std::vector<Eigen::SparseMatrix<double>>& S) {
  Eigen::VectorXd out = Eigen::VectorXd::Zero(chain.size());
  if (S.empty())
    return out;
  Eigen::SparseMatrix<double> total = S.front();
  for (std::size_t k = 1; k < S.size(); ++k)
    total = total + S[k];

  for (int i = 0; i < chain.size(); ++i) {
    const Zai& parent = chain.nuclides()[static_cast<std::size_t>(i)];
    const DecayData* d = chain.decay(parent);
    if (d == nullptr || d->decayConstant == 0.0)
      continue;
    double s = -total.coeff(i, i);
    for (const DecayMode& m : d->modes) {
      if (m.isFission) {
        if (const FissionYields* y = chain.nearestYields(parent, 0.0))
          for (const auto& [product, yield] : y->products)
            s += m.branching * yield * total.coeff(chain.indexOf(product), i);
        continue;
      }
      const std::optional<Zai> daughter = DepletionChain::decayDaughter(parent, m);
      if (!daughter)
        continue;
      const int j = chain.indexOf(*daughter);
      if (j >= 0)
        s += m.branching * total.coeff(j, i);
    }
    out(i) = s;
  }
  return out;
}

namespace {

// dA(row, p) / dsigma for every row one reaction channel of `parent` (at
// matrix index p) touches: -1 for the parent's removal, +1 at a tracked
// target, +yield at each product of the yield table the channel would use.
//
// Accumulated into one map rather than applied term by term so that a product
// which is also the parent nets against the removal instead of replacing it,
// and so the contraction below can be a single walk of one column.
std::unordered_map<int, double> channelGain(const DepletionChain& chain, const Zai& parent, int p,
                                            const ReactionXS& r) {
  std::unordered_map<int, double> gain;
  gain[p] -= 1.0;
  if (r.type == ReactionType::Fission) {
    if (const FissionYields* y = chain.nearestYields(parent, r.energy)) {
      for (const auto& [product, yield] : y->products) {
        const int j = chain.indexOf(product);
        if (j >= 0)
          gain[j] += yield;
      }
    }
  } else if (r.target) {
    const int j = chain.indexOf(*r.target);
    if (j >= 0)
      gain[j] += 1.0;
  }
  return gain;
}

}  // namespace

std::vector<ReactionSensitivity> reactionSensitivities(
    const DepletionSystem& system, const std::vector<Eigen::SparseMatrix<double>>& S,
    const std::vector<double>& flux) {
  if (S.size() != flux.size())
    throw std::invalid_argument("cram: reactionSensitivities needs one flux per S_k");
  constexpr double kBarn = 1e-24;
  const DepletionChain& chain = system.chain();
  const bool wrongSize =
      std::any_of(S.begin(), S.end(), [&chain](const Eigen::SparseMatrix<double>& Sk) {
        return Sk.rows() != chain.size() || Sk.cols() != chain.size();
      });
  if (wrongSize)
    throw std::invalid_argument(
        "cram: reactionSensitivities needs each S_k sized like the system's chain");

  std::vector<ReactionSensitivity> out;
  for (const auto& [parent, rxns] : system.reactions()) {
    const int p = chain.indexOf(parent);
    for (const ReactionXS& r : rxns) {
      double s = 0.0;
      if (p >= 0) {
        const std::unordered_map<int, double> gain = channelGain(chain, parent, p, r);
        for (std::size_t k = 0; k < S.size(); ++k) {
          // Walk column p once and pick out the rows this channel moves.
          // Reading them with coeff() instead would report 0 for an entry S_k
          // never integrated -- indistinguishable from one that integrated to
          // zero -- and, since the removal term is applied either way, would
          // hand back half a derivative with no sign left to trust.
          double dA = 0.0;
          std::size_t found = 0;
          for (Eigen::SparseMatrix<double>::InnerIterator it(S[k], p); it; ++it) {
            const auto g = gain.find(static_cast<int>(it.row()));
            if (g == gain.end())
              continue;
            dA += g->second * it.value();
            ++found;
          }
          if (found != gain.size())
            throw std::invalid_argument(
                "cram: reactionSensitivities: S has no entry for part of the " +
                std::string(reactionName(r.type)) + " channel of " + parent.str() +
                ", so its derivative would be incomplete. S_k carries only the sparsity of "
                "the matrices it was built from, and a channel with no rate on any interval "
                "contributes none; build the interval matrices with intervalMatrices(), which "
                "keeps a structural entry for every channel of the system");
          s += flux[k] * kBarn * dA;
        }
      }
      out.push_back(
          ReactionSensitivity{.parent = parent, .type = r.type, .target = r.target, .dRdSigma = s});
    }
  }
  return out;
}

}  // namespace cram
