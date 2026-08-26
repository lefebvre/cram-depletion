// Example driver for the CRAM depletion library.
//
//   ./deplete                      -> runs a self-contained decay demo (no data files)
//   ./deplete <decay.endf> [t_sec] -> loads ENDF decay data (needs -DWITH_ENDFTK)
//   ./deplete --version            -> prints the library version and exits
//
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "cram/chain.hpp"
#include "cram/cram.hpp"
#include "cram/endf_reader.hpp"
#include "cram/version.hpp"

using namespace cram;

// Build a small illustrative chain by hand and decay it. This exercises the
// matrix assembly + CRAM solver without needing any data files, and the result
// is checked against the analytic two-step Bateman solution.
static void demo() {
  // I-135 (beta-) -> Xe-135 (beta-) -> Cs-135  (half-lives in seconds)
  const Zai I135{53, 135, 0}, Xe135{54, 135, 0}, Cs135{55, 135, 0};
  const double thalf_I = 6.57 * 3600.0;   // 6.57 h
  const double thalf_Xe = 9.14 * 3600.0;  // 9.14 h

  DepletionChain chain;
  chain.add(I135);
  chain.add(Xe135);
  chain.add(Cs135);  // treated as stable here

  chain.setDecay(
      I135, DecayData{.halfLife = thalf_I,
                      .modes = {DecayMode{
                          .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});
  chain.setDecay(
      Xe135, DecayData{.halfLife = thalf_Xe,
                       .modes = {DecayMode{
                           .rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}});

  auto A = chain.decayMatrix();

  Eigen::VectorXd n0(chain.size());
  n0.setZero();
  n0(chain.indexOf(I135)) = 1.0;

  const double t = 12.0 * 3600.0;  // 12 hours
  Eigen::VectorXd n16 = cramSolve(A, n0, t, CramOrder::CRAM16);
  Eigen::VectorXd n48 = cramSolve(A, n0, t, CramOrder::CRAM48);

  // Analytic Bateman for the linear chain a -> b -> c.
  const double la = kLn2 / thalf_I, lb = kLn2 / thalf_Xe;
  const double NA = std::exp(-la * t);
  const double NB = la / (lb - la) * (std::exp(-la * t) - std::exp(-lb * t));
  const double NC = 1.0 - NA - NB;

  std::printf("Decay of 1.0 atom of I-135 over 12 h:\n");
  std::printf("            %-12s %-12s %-12s\n", "I-135", "Xe-135", "Cs-135");
  std::printf("  CRAM16    %-12.6e %-12.6e %-12.6e\n", n16(chain.indexOf(I135)),
              n16(chain.indexOf(Xe135)), n16(chain.indexOf(Cs135)));
  std::printf("  CRAM48    %-12.6e %-12.6e %-12.6e\n", n48(chain.indexOf(I135)),
              n48(chain.indexOf(Xe135)), n48(chain.indexOf(Cs135)));
  std::printf("  analytic  %-12.6e %-12.6e %-12.6e\n", NA, NB, NC);
  std::printf("  mass (CRAM48 sum) = %.15f\n", n48.sum());
}

int main(int argc, char** argv) {
  // Checked before the filename path so `--version` is never mistaken for an
  // ENDF tape. cram::kVersion comes from the generated cram/version.hpp, so
  // this reports the project() version and cannot drift from it.
  if (argc > 1 && (std::strcmp(argv[1], "--version") == 0 || std::strcmp(argv[1], "-v") == 0)) {
    std::printf("deplete (cram) %s\n", cram::kVersion);
    return 0;
  }

  if (argc < 2) {
    demo();
    return 0;
  }

  DepletionChain chain;
  int n = loadDecayData(chain, argv[1]);
  if (n == 0) {
#ifdef WITH_ENDFTK
    std::printf(
        "No MT457 (decay) sections found in '%s'.\n"
        "This driver expects an ENDF/B-VIII decay sublibrary file; NFY (MT454/459)\n"
        "and other tapes are not handled here.\n",
        argv[1]);
#else
    std::printf(
        "This binary was built without WITH_ENDFTK, so ENDF files cannot be read.\n"
        "Rebuild with -DWITH_ENDFTK=ON.\n");
#endif
    return 1;
  }
  std::printf("Loaded decay data for %d nuclides (%d in chain).\n", n, chain.size());

  auto A = chain.decayMatrix();
  Eigen::VectorXd n0(chain.size());
  n0.setZero();
  // Seed each loaded parent with 1.0 atoms so the demo has something to decay.
  // Real use should set the inventory from a problem-specific source.
  for (const auto& z : chain.nuclides()) {
    if (chain.decay(z))
      n0(chain.indexOf(z)) = 1.0;
  }

  const double t = (argc > 2) ? std::stod(argv[2]) : 86400.0;
  Eigen::VectorXd nt = cramSolve(A, n0, t);

  std::printf("Decayed %d-nuclide inventory by %g s.\n", chain.size(), t);
  std::printf("  initial total: %.6e   final total: %.6e\n", n0.sum(), nt.sum());
  std::printf("  nuclide        n0           n(t)         delta\n");
  for (const auto& z : chain.nuclides()) {
    int i = chain.indexOf(z);
    if (std::abs(nt(i)) < 1e-30 && n0(i) == 0.0)
      continue;
    std::printf("  %-12s   %-12.4e %-12.4e %+.3e\n", z.str().c_str(), n0(i), nt(i), nt(i) - n0(i));
  }
  return 0;
}
