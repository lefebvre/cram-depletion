#pragma once
//
// Independent analytic references for the depletion equations, used to check
// the CRAM solver. These are closed-form Bateman solutions, NOT another
// numerical matrix-exponential, so a transcription or algorithm bug in the
// solver cannot hide behind a matching bug in the reference.
//
#include <cmath>
#include <cstddef>
#include <vector>

namespace cram_test {

// Bateman solution for a pure linear decay chain 1 -> 2 -> ... -> n with one
// atom of nuclide 1 at t = 0 and unit branching at every step:
//
//   N_k(t) = (prod_{i=1}^{k-1} lambda_i) *
//            sum_{i=1}^{k} exp(-lambda_i t) / prod_{j=1, j!=i}^{k} (lambda_j - lambda_i)
//
// A stable terminal nuclide is represented with lambda = 0. Requires all
// lambda_i distinct (the closed form has removable singularities otherwise).
inline std::vector<double> batemanLinearChain(const std::vector<double>& lambda, double t) {
  const std::size_t n = lambda.size();
  std::vector<double> N(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) {
    double prodLambda = 1.0;
    for (std::size_t i = 0; i < k; ++i)
      prodLambda *= lambda[i];

    double sum = 0.0;
    for (std::size_t i = 0; i <= k; ++i) {
      double denom = 1.0;
      for (std::size_t j = 0; j <= k; ++j)
        if (j != i)
          denom *= (lambda[j] - lambda[i]);
      sum += std::exp(-lambda[i] * t) / denom;
    }
    N[k] = prodLambda * sum;
  }
  return N;
}

// Two-branch case: parent A (lambda_a) decays to B with branching b_ab and to C
// with branching b_ac; B and C are stable. One atom of A at t = 0.
struct TwoBranchResult {
  double a, b, c;
};
inline TwoBranchResult batemanTwoBranch(double lambda_a, double b_ab, double b_ac, double t) {
  const double eA = std::exp(-lambda_a * t);
  // dB/dt = b_ab*lambda_a*A, integrated with A = exp(-lambda_a t); likewise C.
  return TwoBranchResult{.a = eA, .b = b_ab * (1.0 - eA), .c = b_ac * (1.0 - eA)};
}

}  // namespace cram_test
