#include <benchmark/benchmark.h>

#include "cram/cram.hpp"
#include "cram/cram_solver.hpp"
#include "synthetic_chain.hpp"

namespace {

// Single CRAM solve at varying N. CRAM16 vs CRAM48 are reported separately so
// the cost of the extra poles is visible at a glance.
template <cram::CramOrder Order>
void BM_CramSolve(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const auto chain = cram_bench::buildSyntheticChain(n);
  const double dt = 86400.0;  // one day

  for (auto _ : state) {
    auto out = cram::cramSolve(chain.A, chain.n0, dt, Order);
    benchmark::DoNotOptimize(out);
  }
  state.SetLabel(Order == cram::CramOrder::CRAM16 ? "CRAM16" : "CRAM48");
}
BENCHMARK(BM_CramSolve<cram::CramOrder::CRAM16>)->Arg(256)->Arg(1024)->Arg(1675);
BENCHMARK(BM_CramSolve<cram::CramOrder::CRAM48>)->Arg(256)->Arg(1024)->Arg(1675);

// Multi-step time march using the free function. Each step re-factorizes all
// K poles even though A and dt are unchanged — this is the cost the cached
// solver below replaces.
void BM_DepletionMarch(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const int steps = static_cast<int>(state.range(1));
  const auto chain = cram_bench::buildSyntheticChain(n);
  const double dt = 86400.0;

  for (auto _ : state) {
    Eigen::VectorXd y = chain.n0;
    for (int s = 0; s < steps; ++s)
      y = cram::cramSolve(chain.A, y, dt, cram::CramOrder::CRAM48);
    benchmark::DoNotOptimize(y);
  }
  state.counters["steps"] = steps;
}
BENCHMARK(BM_DepletionMarch)
    ->Args({1675, 20})
    ->Args({1675, 100})
    // ->Args({1675, 500})
    ->Unit(benchmark::kMillisecond);

// Same march, but uses CramSolver: prepare() once, then N apply() calls. The
// per-iteration cost should drop to just K complex sparse solves (24 for
// CRAM48), with the factorize cost amortized into the first iteration.
void BM_CramSolverCached(benchmark::State& state) {
  const int n = static_cast<int>(state.range(0));
  const int steps = static_cast<int>(state.range(1));
  const auto chain = cram_bench::buildSyntheticChain(n);
  const double dt = 86400.0;

  for (auto _ : state) {
    cram::CramSolver solver(cram::CramOrder::CRAM48);
    solver.prepare(chain.A, dt);
    Eigen::VectorXd y = chain.n0;
    for (int s = 0; s < steps; ++s)
      y = solver.apply(y);
    benchmark::DoNotOptimize(y);
  }
  state.counters["steps"] = steps;
}
BENCHMARK(BM_CramSolverCached)
    ->Args({1675, 20})
    ->Args({1675, 100})
    // ->Args({1675, 500})
    ->Unit(benchmark::kMillisecond);

}  // namespace
