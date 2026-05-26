#include <benchmark/benchmark.h>

#include "cram/cram.hpp"
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

// Multi-step time march. Calls cramSolve repeatedly with the same A, mirroring
// a real depletion run. Each step is an independent matrix exponential so
// nothing carries across calls today — this benchmark is the place to track
// improvements like caching the symbolic factorization across steps.
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
    ->Args({1675, 5})
    ->Args({1675, 20})
    ->Unit(benchmark::kMillisecond);

}  // namespace
