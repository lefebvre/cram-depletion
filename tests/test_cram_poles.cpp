#include "cram/cram_poles.hpp"

// The header's promise is that code which cannot compile Eigen can include it.
// Nothing else in this file brings Eigen in, so an Eigen macro seen here came
// from cram_poles.hpp.
#if defined(EIGEN_WORLD_VERSION) || defined(EIGEN_MAJOR_VERSION)
#error "cram/cram_poles.hpp must not include Eigen"
#endif

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>

using cram::CramOrder;
using cram::CramPoles;
using cram::cramPoles;

namespace {

// cramPoles() is constexpr, so a caller can size or consume a table at compile
// time.
static_assert(cramPoles(CramOrder::CRAM16).theta.size() == 8);
static_assert(cramPoles(CramOrder::CRAM48).theta.size() == 24);

// The scalar IPF form evaluated straight from a table: exp(z) for the 1x1
// "matrix" z = lambda t.
double ipfScalar(const CramPoles& p, double z) {
  double y = 1.0;
  for (std::size_t l = 0; l < p.theta.size(); ++l) {
    const std::complex<double> x = y / (z - p.theta[l]);
    y += 2.0 * (p.alpha[l] * x).real();
  }
  return y * p.alpha0;
}

}  // namespace

TEST(CramPoles, EachOrderHasHalfItsDegreeInPoles) {
  for (CramOrder order : {CramOrder::CRAM16, CramOrder::CRAM48}) {
    const CramPoles p = cramPoles(order);
    EXPECT_EQ(p.theta.size(), static_cast<std::size_t>(order) / 2);
    EXPECT_EQ(p.alpha.size(), p.theta.size());
    EXPECT_GT(p.alpha0, 0.0);
  }
}

// Documented fallback, shared with cramSolve() and CramSolver.
TEST(CramPoles, AnyOtherOrderReadsAsCram48) {
  const CramPoles fallback = cramPoles(static_cast<CramOrder>(0));
  const CramPoles cram48 = cramPoles(CramOrder::CRAM48);
  EXPECT_EQ(fallback.theta.data(), cram48.theta.data());
  EXPECT_EQ(fallback.alpha.data(), cram48.alpha.data());
  EXPECT_EQ(fallback.alpha0, cram48.alpha0);
}

// |lambda t - theta| >= Im(theta) for any real lambda t, so this bounds the
// divisor of a substitution on a real triangular matrix away from zero. The
// smallest imaginary part is 1.194 in both orders.
TEST(CramPoles, EveryPoleIsOffTheRealAxis) {
  for (CramOrder order : {CramOrder::CRAM16, CramOrder::CRAM48}) {
    for (const std::complex<double>& theta : cramPoles(order).theta)
      EXPECT_GT(theta.imag(), 1.0) << "order " << static_cast<int>(order) << " theta " << theta;
  }
}

// The worst error measured on a geometric sweep of [-1e4, 0] is 1.3e-15 for
// CRAM16 and 2.4e-15 for CRAM48; the tolerance leaves about 4x headroom.
TEST(CramPoles, ScalarFormReproducesExp) {
  constexpr double kTol = 1e-14;
  for (CramOrder order : {CramOrder::CRAM16, CramOrder::CRAM48}) {
    const CramPoles p = cramPoles(order);
    for (double z : {0.0, -1e-8, -0.03, -1.0, -10.0, -100.0, -1e3, -1e4})
      EXPECT_NEAR(ipfScalar(p, z), std::exp(z), kTol)
          << "order " << static_cast<int>(order) << " z " << z;
  }
}
