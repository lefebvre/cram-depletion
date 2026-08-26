#pragma once
//
// A small but representative thermal-pin depletion fixture used by the
// integrator, depletion and adjoint tests. It is NOT the real VERA pin (that
// needs transport-generated cross sections, see validation/openmc) -- it is a
// hand-built chain with plausible one-group cross sections that exercises the
// same machinery: fissile depletion, Pu build-in, an I-135 -> Xe-135 -> Cs-135
// fission-product chain with a strong Xe-135 absorber, all under constant-power
// flux normalization so the burnup matrix is genuinely composition-dependent.
//
#include <Eigen/Core>
#include <optional>
#include <vector>

#include "cram/chain.hpp"
#include "cram/deplete.hpp"
#include "cram/nuclide.hpp"
#include "cram/reaction.hpp"

namespace cram_test {

inline constexpr cram::Zai kU235{.z = 92, .a = 235, .i = 0};
inline constexpr cram::Zai kU236{.z = 92, .a = 236, .i = 0};
inline constexpr cram::Zai kU238{.z = 92, .a = 238, .i = 0};
inline constexpr cram::Zai kU239{.z = 92, .a = 239, .i = 0};
inline constexpr cram::Zai kNp239{.z = 93, .a = 239, .i = 0};
inline constexpr cram::Zai kPu239{.z = 94, .a = 239, .i = 0};
inline constexpr cram::Zai kPu240{.z = 94, .a = 240, .i = 0};
inline constexpr cram::Zai kPu241{.z = 94, .a = 241, .i = 0};
inline constexpr cram::Zai kPu242{.z = 94, .a = 242, .i = 0};
inline constexpr cram::Zai kAm241{.z = 95, .a = 241, .i = 0};
inline constexpr cram::Zai kI135{.z = 53, .a = 135, .i = 0};
inline constexpr cram::Zai kXe135{.z = 54, .a = 135, .i = 0};
inline constexpr cram::Zai kXe136{.z = 54, .a = 136, .i = 0};
inline constexpr cram::Zai kCs135{.z = 55, .a = 135, .i = 0};

constexpr double kFissionQ = 200.0e6;  // eV per fission (~200 MeV)

// Build the fixture chain: register nuclides, decay data, and fission yields.
inline cram::DepletionChain buildPinChain() {
  using namespace cram;
  DepletionChain chain;
  for (const Zai& z : {kU235, kU236, kU238, kU239, kNp239, kPu239, kPu240, kPu241, kPu242, kAm241,
                       kI135, kXe135, kXe136, kCs135})
    chain.add(z);

  const auto betaMinus = [](double thalf) {
    return DecayData{
        .halfLife = thalf,
        .modes = {DecayMode{.rtyp = 1.0, .branching = 1.0, .finalState = 0, .isFission = false}}};
  };
  chain.setDecay(kU239, betaMinus(1407.0));     // 23.45 min -> Np239
  chain.setDecay(kNp239, betaMinus(203558.0));  // 2.356 d  -> Pu239
  chain.setDecay(kI135, betaMinus(23652.0));    // 6.57 h   -> Xe135
  chain.setDecay(kXe135, betaMinus(32904.0));   // 9.14 h   -> Cs135
  chain.setDecay(kPu241, betaMinus(4.52e8));    // 14.3 y   -> Am241

  // Thermal independent fission yields (approximate, direct-to-chain members).
  chain.addFissionYields(
      kU235, FissionYields{.energy = 0.0253, .products = {{kI135, 0.0292}, {kXe135, 0.00254}}});
  chain.addFissionYields(
      kPu239, FissionYields{.energy = 0.0253, .products = {{kI135, 0.0604}, {kXe135, 0.0105}}});
  return chain;
}

inline cram::ReactionXS fission(double sigma) {
  return cram::ReactionXS{
      .type = cram::ReactionType::Fission, .target = std::nullopt, .sigma = sigma, .q = kFissionQ};
}

inline cram::ReactionXS capture(const cram::Zai& product, double sigma) {
  return cram::ReactionXS{.type = cram::ReactionType::NGamma, .target = product, .sigma = sigma};
}

// Configure the one-group cross sections (barn) on a system built from
// buildPinChain(). Representative thermal-pin values.
inline void configurePinReactions(cram::DepletionSystem& sys) {
  sys.setReactions(kU235, {fission(38.0), capture(kU236, 9.0)});
  sys.setReactions(kU238, {fission(0.05), capture(kU239, 0.9)});
  sys.setReactions(kPu239, {fission(105.0), capture(kPu240, 59.0)});
  sys.setReactions(kPu240, {capture(kPu241, 290.0)});
  sys.setReactions(kPu241, {fission(130.0), capture(kPu242, 41.0)});
  sys.setReactions(kXe135, {capture(kXe136, 2.6e6)});
}

// Initial number densities (atom/b-cm) for ~3.1 w/o UO2.
inline Eigen::VectorXd initialPinComposition(const cram::DepletionChain& chain) {
  Eigen::VectorXd n0(chain.size());
  n0.setZero();
  n0(chain.indexOf(kU235)) = 7.0e-4;
  n0(chain.indexOf(kU238)) = 2.2e-2;
  return n0;
}

}  // namespace cram_test
