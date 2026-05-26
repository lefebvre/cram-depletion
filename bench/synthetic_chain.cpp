#include "synthetic_chain.hpp"

#include <cmath>
#include <random>
#include <vector>

namespace cram_bench {

SyntheticChain buildSyntheticChain(int n, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> logLambda(-9.0, -1.0);  // 1e-9 .. 1e-1
  std::uniform_int_distribution<int> nModes(1, 3);
  std::uniform_int_distribution<int> hop(1, 4);
  std::uniform_real_distribution<double> branch(0.1, 1.0);

  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n) * 4);

  for (int i = 0; i < n - 1; ++i) {
    const double lambda = std::pow(10.0, logLambda(rng));
    const int modes = nModes(rng);

    std::vector<double> br(modes);
    double sum = 0.0;
    for (int k = 0; k < modes; ++k) {
      br[k] = branch(rng);
      sum += br[k];
    }
    for (double& b : br)
      b /= sum;

    t.emplace_back(i, i, -lambda);
    for (int k = 0; k < modes; ++k) {
      int j = i + hop(rng);
      if (j >= n)
        j = n - 1;
      t.emplace_back(j, i, lambda * br[k]);
    }
  }

  SyntheticChain c;
  c.A.resize(n, n);
  c.A.setFromTriplets(t.begin(), t.end());  // sums duplicates
  c.A.makeCompressed();

  c.n0 = Eigen::VectorXd::Zero(n);
  c.n0(0) = 1.0;
  return c;
}

}  // namespace cram_bench
