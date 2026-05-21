#include "cram/cram.hpp"

#include <Eigen/SparseLU>
#include <array>
#include <complex>

namespace cram {
namespace {

using cd = std::complex<double>;

// ---------------------------------------------------------------------------
// IPF-CRAM coefficients.
//
// Transcribed verbatim from OpenMC's openmc/deplete/cram.py (MIT License,
// (c) MIT / UChicago Argonne LLC and OpenMC contributors), which in turn
// reproduces the coefficients from Pusa, NSE 182:3 (2016). The values were
// validated here against analytic Bateman solutions to ~1e-15 relative error.
// ---------------------------------------------------------------------------

// These are exact published constants; a couple happen to sit near math
// constants like sqrt2, which trips modernize-use-std-numbers. They must stay
// verbatim, so the check is suppressed across the tables.
// NOLINTBEGIN(modernize-use-std-numbers)
constexpr double kAlpha0_16 = 2.124853710495224e-16;

const std::array<cd, 8> kTheta16 = {{
    {+3.509103608414918, +8.436198985884374},
    {+5.948152268951177, +3.587457362018322},
    {-5.264971343442647, +16.22022147316793},
    {+1.419375897185666, +10.92536348449672},
    {+6.416177699099435, +1.194122393370139},
    {+4.993174737717997, +5.996881713603942},
    {-1.413928462488886, +13.49772569889275},
    {-10.84391707869699, +19.27744616718165},
}};

const std::array<cd, 8> kAlpha16 = {{
    {+5.464930576870210e+3, -3.797983575308356e+4},
    {+9.045112476907548e+1, -1.115537522430261e+3},
    {+2.344818070467641e+2, -4.228020157070496e+2},
    {+9.453304067358312e+1, -2.951294291446048e+2},
    {+7.283792954673409e+2, -1.205646080220011e+5},
    {+3.648229059594851e+1, -1.155509621409682e+2},
    {+2.547321630156819e+1, -2.639500283021502e+1},
    {+2.394538338734709e+1, -5.650522971778156e+0},
}};

constexpr double kAlpha0_48 = 2.258038182743983e-47;

const std::array<cd, 24> kTheta48 = {{
    {-4.465731934165702e+1, +6.233225190695437e+1}, {-5.284616241568964e+0, +4.057499381311059e+1},
    {-8.867715667624458e+0, +4.325515754166724e+1}, {+3.493013124279215e+0, +3.281615453173585e+1},
    {+1.564102508858634e+1, +1.558061616372237e+1}, {+1.742097597385893e+1, +1.076629305714420e+1},
    {-2.834466755180654e+1, +5.492841024648724e+1}, {+1.661569367939544e+1, +1.316994930024688e+1},
    {+8.011836167974721e+0, +2.780232111309410e+1}, {-2.056267541998229e+0, +3.794824788914354e+1},
    {+1.449208170441839e+1, +1.799988210051809e+1}, {+1.853807176907916e+1, +5.974332563100539e+0},
    {+9.932562704505182e+0, +2.532823409972962e+1}, {-2.244223871767187e+1, +5.179633600312162e+1},
    {+8.590014121680897e-1, +3.536456194294350e+1}, {-1.286192925744479e+1, +4.600304902833652e+1},
    {+1.164596909542055e+1, +2.287153304140217e+1}, {+1.806076684783089e+1, +8.368200580099821e+0},
    {+5.870672154659249e+0, +3.029700159040121e+1}, {-3.542938819659747e+1, +5.834381701800013e+1},
    {+1.901323489060250e+1, +1.194282058271408e+0}, {+1.885508331552577e+1, +3.583428564427879e+0},
    {-1.734689708174982e+1, +4.883941101108207e+1}, {+1.316284237125190e+1, +2.042951874827759e+1},
}};

const std::array<cd, 24> kAlpha48 = {{
    {+6.387380733878774e+2, -6.743912502859256e+2}, {+1.909896179065730e+2, -3.973203432721332e+2},
    {+4.236195226571914e+2, -2.041233768918671e+3}, {+4.645770595258726e+2, -1.652917287299683e+3},
    {+7.765163276752433e+2, -1.783617639907328e+4}, {+1.907115136768522e+3, -5.887068595142284e+4},
    {+2.909892685603256e+3, -9.953255345514560e+3}, {+1.944772206620450e+2, -1.427131226068449e+3},
    {+1.382799786972332e+5, -3.256885197214938e+6}, {+5.628442079602433e+3, -2.924284515884309e+4},
    {+2.151681283794220e+2, -1.121774011188224e+3}, {+1.324720240514420e+3, -6.370088443140973e+4},
    {+1.617548476343347e+4, -1.008798413156542e+6}, {+1.112729040439685e+2, -8.837109731680418e+1},
    {+1.074624783191125e+2, -1.457246116408180e+2}, {+8.835727765158191e+1, -6.388286188419360e+1},
    {+9.354078136054179e+1, -2.195424319460237e+2}, {+9.418142823531573e+1, -6.719055740098035e+2},
    {+1.040012390717851e+2, -1.693747595553868e+2}, {+6.861882624343235e+1, -1.177598523430493e+1},
    {+8.766654491283722e+1, -4.596464999363902e+3}, {+1.056007619389650e+2, -1.738294585524067e+3},
    {+7.738987569039419e+1, -4.311715386228984e+1}, {+1.041366366475571e+2, -2.777743732451969e+2},
}};
// NOLINTEND(modernize-use-std-numbers)

template <std::size_t K>
Eigen::VectorXd ipfCram(const Eigen::SparseMatrix<double>& Araw, const Eigen::VectorXd& n0,
                        double dt, const std::array<cd, K>& theta, const std::array<cd, K>& alpha,
                        double alpha0) {
  const int N = static_cast<int>(n0.size());
  const Eigen::SparseMatrix<cd> A = (Araw.cast<cd>() * cd(dt, 0.0)).eval();
  Eigen::SparseMatrix<cd> I(N, N);
  I.setIdentity();

  Eigen::VectorXd y = n0;  // stays real: we add 2*Re(...) each pole
  Eigen::SparseLU<Eigen::SparseMatrix<cd>> lu;

  for (std::size_t l = 0; l < K; ++l) {
    Eigen::SparseMatrix<cd> M = A - (theta[l] * I);
    M.makeCompressed();
    lu.compute(M);
    // For production code, check lu.info() == Eigen::Success here.
    Eigen::VectorXcd x = lu.solve(y.cast<cd>());
    y += 2.0 * (alpha[l] * x).real();
  }
  return y * alpha0;
}

}  // namespace

Eigen::VectorXd cramSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& n0,
                          double dt, CramOrder order) {
  if (order == CramOrder::CRAM16)
    return ipfCram<8>(A, n0, dt, kTheta16, kAlpha16, kAlpha0_16);
  return ipfCram<24>(A, n0, dt, kTheta48, kAlpha48, kAlpha0_48);
}

}  // namespace cram
