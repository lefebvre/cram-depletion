// Validation test reproducing a Pusa-2016 actinide decay case
// (Pu-241 -> Am-241 -> Np-237 -> Pa-233 -> U-233) against real
// ENDF/B-VIII.0 decay tapes fetched at configure time by
// cmake/FetchEndfActinides.cmake.
//
// Reference: M. Pusa, "Higher-Order Chebyshev Rational Approximation
// Method and Application to Burnup Equations," NSE 182:3, 297-318 (2016).
//
// Current assertions use robust kinetic invariants (half-life behaviour,
// 10-half-life tail, mass non-creation). The TODO at the bottom marks
// where exact reference inventories from the paper should be inserted
// to upgrade this from verification-grade to validation-grade.
//
// Built only with WITH_ENDFTK.
#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/endf_reader.hpp"

using namespace cram;

namespace {

struct Step {
  Zai zai;
  const char* file;
  double halfLifeReference;  // seconds, for a sanity check on the loaded tape
};

constexpr double kYear = 365.25 * 86400.0;
constexpr double kDay = 86400.0;

// Half-lives are from standard references; ENDF/B-VIII.0 may disagree in
// the last decimal, so the test tolerates 0.1% drift.
const std::array<Step, 5> kChain = {{
    {{.z = 94, .a = 241, .i = 0}, "decay_3561_94-Pu-241.dat", 14.290 * kYear},
    {{.z = 95, .a = 241, .i = 0}, "decay_3578_95-Am-241.dat", 432.6 * kYear},
    {{.z = 93, .a = 237, .i = 0}, "decay_3537_93-Np-237.dat", 2.144e6 * kYear},
    {{.z = 91, .a = 233, .i = 0}, "decay_3489_91-Pa-233.dat", 26.975 * kDay},
    {{.z = 92, .a = 233, .i = 0}, "decay_3513_92-U-233.dat", 1.591e5 * kYear},
}};

}  // namespace

TEST(PusaActinideValidation, Pu241ActinideChain) {
  DepletionChain chain;
  for (const auto& step : kChain) {
    int n = loadDecayData(chain, step.file);
    ASSERT_EQ(n, 1) << "expected one MT457 in " << step.file;
    const DecayData* d = chain.decay(step.zai);
    ASSERT_NE(d, nullptr) << step.zai.str();
    EXPECT_NEAR(d->halfLife, step.halfLifeReference, step.halfLifeReference * 1e-3)
        << "half-life drift for " << step.zai.str();
  }

  auto A = chain.decayMatrix();
  Eigen::VectorXd n0(chain.size());
  n0.setZero();
  const int iPu241 = chain.indexOf(kChain[0].zai);
  ASSERT_GE(iPu241, 0);
  n0(iPu241) = 1.0;

  const double tHalfPu241 = chain.decay(kChain[0].zai)->halfLife;

  // At t = T1/2(Pu-241), Pu-241 fraction is exactly 0.5.
  {
    Eigen::VectorXd nt = cramSolve(A, n0, tHalfPu241, CramOrder::CRAM48);
    EXPECT_NEAR(nt(iPu241), 0.5, 1e-9) << "Pu-241 at one half-life";
  }

  // At t = 10 * T1/2(Pu-241), Pu-241 fraction is 2^-10 ~= 9.77e-4.
  {
    Eigen::VectorXd nt = cramSolve(A, n0, 10.0 * tHalfPu241, CramOrder::CRAM48);
    EXPECT_NEAR(nt(iPu241), std::pow(2.0, -10), 1e-9);
    EXPECT_GE(nt(iPu241), 0.0) << "negative inventory";
  }

  // Mass conservation: the only path out of the modelled chain is the
  // U-233 -> Th-229 alpha decay, but Th-229 is added by loadDecayData as a
  // zero-decay sink, so total atom count tracks the initial 1.0. ENDF/B-VIII.0
  // branching ratios are published with finite precision and do not always
  // sum to exactly 1.0 per parent (a 1e-6 excess is typical), so the bound is
  // set accordingly rather than at machine epsilon.
  {
    Eigen::VectorXd nt = cramSolve(A, n0, 1.0e3 * tHalfPu241, CramOrder::CRAM48);
    EXPECT_NEAR(nt.sum(), 1.0, 1e-5) << "mass non-conservation beyond ENDF BR precision";
  }

  // TODO(validation): replace the invariants above with exact reference
  // inventories tabulated in Pusa, NSE 182:3 (2016) (specific table to be
  // identified). Until those values are inserted here this test is
  // verification-grade, not full validation-grade.
}
