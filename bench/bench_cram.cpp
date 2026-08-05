// Benchmarks over both matrix regimes this library serves.
//
//   decay-only  acyclic; every transition runs forward. What
//               DepletionChain::decayMatrix() builds, and what the ENDF MT457
//               path and cram-apps/deplete.cpp operate on. SparseLU finds an
//               elimination order with little fill-in, so these are the cheap
//               ones.
//   burnup      cyclic; (n,gamma) walks the chain up while decay walks it back
//               down, plus fission products linking high indices to low ones.
//               Fill-in dominates the factorization.
//
// Both are measured because they cost substantially different amounts and an
// optimization that helps one need not help the other. Reporting only one of
// them is how the suite previously came to describe an acyclic matrix as
// "representative of a real depletion chain".
#include <benchmark/benchmark.h>

#include "acyclic_chain.hpp"
#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/cram_solver.hpp"
#include "cyclic_chain.hpp"
#include "synthetic_depletion_chain.hpp"

namespace {

constexpr double kDay = 86400.0;

enum class Shape { DecayOnly, Burnup };

template <Shape S>
auto buildChain(int n) {
  if constexpr (S == Shape::DecayOnly) {
    return cram_test::buildAcyclicChain(n);
  } else {
    return cram_test::buildCyclicChain(n);
  }
}

constexpr const char* shapeName(Shape s) {
  return s == Shape::DecayOnly ? "decay-only" : "burnup";
}

// One-shot solve: factorize all K poles and solve, the cost paid per call by
// cramSolve(). CRAM16 and CRAM48 are reported separately so the price of the
// extra poles is visible.
template <Shape S, cram::CramOrder Order>
void BM_CramSolve(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto chain = buildChain<S>(n);

  for (auto _ : state) {
    auto out = cram::cramSolve(chain.A, chain.n0, kDay, Order);
    benchmark::DoNotOptimize(out);
  }
  state.SetLabel(std::string(shapeName(S)) +
                 (Order == cram::CramOrder::CRAM16 ? "/CRAM16" : "/CRAM48"));
}
BENCHMARK(BM_CramSolve<Shape::DecayOnly, cram::CramOrder::CRAM16>)->Arg(256)->Arg(1024)->Arg(1675);
BENCHMARK(BM_CramSolve<Shape::DecayOnly, cram::CramOrder::CRAM48>)->Arg(256)->Arg(1024)->Arg(1675);
BENCHMARK(BM_CramSolve<Shape::Burnup, cram::CramOrder::CRAM16>)->Arg(256)->Arg(1024)->Arg(1675);
BENCHMARK(BM_CramSolve<Shape::Burnup, cram::CramOrder::CRAM48>)->Arg(256)->Arg(1024)->Arg(1675);

// The two halves of the cached solver, measured apart. Reporting prepare+march
// as one number hides which half a change moved, and the halves have opposite
// parallelism: the K factorizations are independent of one another, while
// apply() is strictly sequential because the IPF form feeds each pole's result
// into the next.
template <Shape S>
void BM_CramSolverPrepare(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto chain = buildChain<S>(n);

  for (auto _ : state) {
    cram::CramSolver solver(cram::CramOrder::CRAM48);
    solver.prepare(chain.A, kDay);
    benchmark::DoNotOptimize(solver);
  }
  state.SetLabel(shapeName(S));
}
BENCHMARK(BM_CramSolverPrepare<Shape::DecayOnly>)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(1675)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_CramSolverPrepare<Shape::Burnup>)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(1675)
    ->Unit(benchmark::kMillisecond);

// The many-region sweep: one solver reused across regions that share a chain
// topology, so every prepare() after the first reuses the symbolic analysis and
// pays only the numeric factorization. Compare against BM_CramSolverPrepare
// above, which constructs a fresh solver each iteration and therefore always
// takes the cold path. The gap is the symbolic analysis.
template <Shape S>
void BM_CramSolverRepreparePreanalyzed(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto chain = buildChain<S>(n);
  cram::CramSolver solver(cram::CramOrder::CRAM48);
  solver.prepare(chain.A, kDay);  // cold prepare kept out of the timed region

  for (auto _ : state) {
    solver.prepare(chain.A, kDay);
    benchmark::DoNotOptimize(solver);
  }
  if (!solver.reusedSymbolicAnalysis())
    state.SkipWithError("expected the symbolic analysis to be reused");
  state.SetLabel(shapeName(S));
}
BENCHMARK(BM_CramSolverRepreparePreanalyzed<Shape::DecayOnly>)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(1675)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_CramSolverRepreparePreanalyzed<Shape::Burnup>)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(1675)
    ->Unit(benchmark::kMillisecond);

// Solves from the same input vector every iteration rather than marching, so
// the inventory cannot decay into denormals and skew the timing.
template <Shape S>
void BM_CramSolverApply(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto chain = buildChain<S>(n);
  cram::CramSolver solver(cram::CramOrder::CRAM48);
  solver.prepare(chain.A, kDay);

  for (auto _ : state) {
    Eigen::VectorXd out = solver.apply(chain.n0);
    benchmark::DoNotOptimize(out);
  }
  state.SetLabel(shapeName(S));
}
BENCHMARK(BM_CramSolverApply<Shape::DecayOnly>)->Arg(256)->Arg(1024)->Arg(1675);
BENCHMARK(BM_CramSolverApply<Shape::Burnup>)->Arg(256)->Arg(1024)->Arg(1675);

// End-to-end march via the free function: every step refactorizes all K poles
// even though A and dt never change. This is the cost CramSolver removes.
template <Shape S>
void BM_DepletionMarch(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const int steps = static_cast<int>(state.range(1));
  const auto chain = buildChain<S>(n);

  for (auto _ : state) {
    Eigen::VectorXd y = chain.n0;
    for (int s = 0; s < steps; ++s)
      y = cram::cramSolve(chain.A, y, kDay, cram::CramOrder::CRAM48);
    benchmark::DoNotOptimize(y);
  }
  state.counters["steps"] = steps;
  state.SetLabel(shapeName(S));
}
BENCHMARK(BM_DepletionMarch<Shape::DecayOnly>)
    ->Args({1675, 20})
    ->Args({1675, 100})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_DepletionMarch<Shape::Burnup>)
    ->Args({1675, 20})
    ->Args({1675, 100})
    ->Unit(benchmark::kMillisecond);

// The same march through CramSolver: prepare() once, then N apply() calls.
// Decay-only marches over many sub-steps are the case CramSolver's header calls
// out as its reason to exist, so that regime is measured too.
template <Shape S>
void BM_CramSolverCached(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const int steps = static_cast<int>(state.range(1));
  const auto chain = buildChain<S>(n);

  for (auto _ : state) {
    cram::CramSolver solver(cram::CramOrder::CRAM48);
    solver.prepare(chain.A, kDay);
    Eigen::VectorXd y = chain.n0;
    for (int s = 0; s < steps; ++s)
      y = solver.apply(y);
    benchmark::DoNotOptimize(y);
  }
  state.counters["steps"] = steps;
  state.SetLabel(shapeName(S));
}
BENCHMARK(BM_CramSolverCached<Shape::DecayOnly>)
    ->Args({1675, 20})
    ->Args({1675, 100})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_CramSolverCached<Shape::Burnup>)
    ->Args({1675, 20})
    ->Args({1675, 100})
    ->Unit(benchmark::kMillisecond);

// --- matrix assembly -------------------------------------------------------
// Assembly runs as often as the solver does: once per depletion region, again
// every burnup step. These cover the two paths -- decay transitions, which walk
// every mode of every nuclide, and the fission source, which walks a yield
// table of ~1000 products per parent.

void BM_DecayMatrixAssembly(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto fixture = cram_test::buildSyntheticDepletionChain(n);

  for (auto _ : state) {
    auto A = fixture.chain.decayMatrix();
    benchmark::DoNotOptimize(A);
  }
  state.counters["nuclides"] = fixture.chain.size();
}
BENCHMARK(BM_DecayMatrixAssembly)->Arg(256)->Arg(1024)->Arg(1675);

void BM_FissionSourceAssembly(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto fixture = cram_test::buildSyntheticDepletionChain(n);

  for (auto _ : state) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (const auto& parent : fixture.fissionParents)
      fixture.chain.addFissionSource(triplets, parent, 1.0e-8, 0.0253);
    benchmark::DoNotOptimize(triplets);
  }
  state.counters["parents"] = static_cast<double>(fixture.fissionParents.size());
}
BENCHMARK(BM_FissionSourceAssembly)->Arg(256)->Arg(1024)->Arg(1675);

}  // namespace
