#pragma once
//
// Nuclide identity (Z, A, isomeric state) and the ENDF decay-mode
// transition rules needed to figure out the daughter of a decay.
//
#include <cstdint>
#include <format>
#include <functional>
#include <string>
#include <vector>

namespace cram {

// A nuclide is identified by (Z, A, I) where I is the isomeric/metastable
// level (0 = ground state, 1 = first metastable, ...). We pack it into a
// single integer key for fast hashing:  ZAI = Z*10000 + A*10 + I.
struct Zai {
  int z = 0;  // proton number
  int a = 0;  // mass number
  int i = 0;  // isomeric state (LISO/RFS)

  constexpr std::int64_t key() const {
    return (std::int64_t(z) * 10000) + (std::int64_t(a) * 10) + i;
  }
  constexpr bool operator==(const Zai& o) const { return key() == o.key(); }

  std::string str() const {
    std::string s = std::to_string(z) + "-" + std::to_string(a);
    if (i != 0)
      s += "m" + std::to_string(i);
    return s;
  }
};

struct ZaiHash {
  std::size_t operator()(const Zai& z) const noexcept { return std::hash<std::int64_t>{}(z.key()); }
};

// ENDF MT457 RTYP single-digit decay codes.
enum class DecayParticle {
  Gamma = 0,            // 0  isomeric transition photon (no Z/A change)
  BetaMinus = 1,        // 1  beta-           Z+1
  ElectronCapture = 2,  // 2  beta+ / EC      Z-1
  IsomericTransition = 3,
  Alpha = 4,    // 4  alpha           Z-2, A-4
  Neutron = 5,  // 5  neutron         A-1
  SpontaneousFission = 6,
  Proton = 7,  // 7  proton          Z-1, A-1
};

// ENDF stores a (possibly multi-step) decay mode in the single real number
// RTYP: the integer part is the first emitted particle and each digit after
// the decimal point is the next particle, in order. e.g.
//   1.0  -> beta-
//   1.5  -> beta- followed by delayed neutron
//   1.4  -> beta- followed by alpha
//   2.7  -> EC/beta+ followed by proton
// Returns the ordered list of single-digit codes.
inline std::vector<int> decayParticleSequence(double rtyp) {
  std::string s = std::format("{:.8g}", rtyp);
  std::vector<int> seq;
  auto dot = s.find('.');
  std::string ip = (dot == std::string::npos) ? s : s.substr(0, dot);
  if (!ip.empty())
    seq.push_back(std::stoi(ip));
  if (dot != std::string::npos)
    for (char c : s.substr(dot + 1))
      if (c >= '0' && c <= '9')
        seq.push_back(c - '0');
  return seq;
}

// Apply one RTYP chain to a parent nuclide and return the daughter.
// `finalState` is the ENDF RFS (isomeric state of the *final* daughter).
// If the chain contains spontaneous fission, `isFission` is set true and the
// returned Zai is meaningless (fission products come from the SFY/NFY tables).
inline Zai applyDecay(const Zai& parent, double rtyp, int finalState, bool& isFission) {
  Zai d = parent;
  d.i = 0;
  isFission = false;
  for (int code : decayParticleSequence(rtyp)) {
    switch (static_cast<DecayParticle>(code)) {
      case DecayParticle::Gamma:
      case DecayParticle::IsomericTransition:
        break;
      case DecayParticle::BetaMinus:
        d.z += 1;
        break;
      case DecayParticle::ElectronCapture:
        d.z -= 1;
        break;
      case DecayParticle::Alpha:
        d.z -= 2;
        d.a -= 4;
        break;
      case DecayParticle::Neutron:
        d.a -= 1;
        break;
      case DecayParticle::Proton:
        d.z -= 1;
        d.a -= 1;
        break;
      case DecayParticle::SpontaneousFission:
        isFission = true;
        break;
      default:
        break;  // unknown code: leave unchanged
    }
  }
  d.i = finalState;
  return d;
}

}  // namespace cram
