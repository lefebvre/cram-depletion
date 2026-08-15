#pragma once
//
// Nuclide identity (Z, A, isomeric state) and the ENDF decay-mode
// transition rules needed to figure out the daughter of a decay.
//
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>

namespace cram {

// A nuclide is identified by (Z, A, I) where I is the isomeric/metastable
// level (0 = ground state, 1 = first metastable, ...). We pack it into a
// single integer key for fast hashing:  ZAI = Z*10000 + A*10 + I.
//
// Z and A deliberately have no default member initializer. No nuclide has Z = 0
// or A = 0, so zero is not an "unset" nuclide but a phantom one: it packs to key
// 0, which add() registers and indexes like any other, and a forgotten field
// therefore yields a chain entry that looks structurally fine and matches
// nothing real. Leaving the initializers out makes -Wmissing-field-initializers
// demand both at each braced initialization. `i` keeps its initializer because
// the ground state is genuinely the right default for an unstated isomeric
// level, which is what lets the common Zai{92, 235} stay short. Construct with
// braces; a bare `Zai z;` leaves Z and A indeterminate and is caught only by
// -Wuninitialized, which cannot see every path.
struct Zai {
  int z;      // proton number
  int a;      // mass number
  int i = 0;  // isomeric state (LISO/RFS); 0 = ground state

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

// Provided so callers can key their own std::unordered_map/set on Zai directly.
// DepletionChain deliberately does not use it -- it keys on the packed int64
// instead -- so this is API for consumers, not an internal detail, and is
// exercised only by its own test.
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
//
// Fixed capacity rather than a std::vector: this is called once per decay mode
// every time a burnup matrix is assembled, and a model with many depletion
// regions assembles one matrix per region. A full ENDF/B-VIII decay sublibrary
// runs to ~3800 nuclides with a few modes each, so a heap allocation here costs
// on the order of 10^4 allocations per matrix build, repeated per region. The
// sequence can never exceed nine entries -- the integer part plus the eight
// fractional digits kScale resolves -- so it fits in the object.
class DecaySequence {
public:
  void push(int code) { codes_[static_cast<std::size_t>(count_++)] = code; }
  const int* begin() const { return codes_.data(); }
  const int* end() const { return codes_.data() + count_; }
  int size() const { return count_; }
  bool empty() const { return count_ == 0; }
  int operator[](int i) const { return codes_[static_cast<std::size_t>(i)]; }

private:
  static constexpr int kMaxCodes = 9;
  std::array<int, kMaxCodes> codes_{};
  int count_ = 0;
};

inline DecaySequence decayParticleSequence(double rtyp) {
  DecaySequence seq;
  if (!(rtyp >= 0.0))
    return seq;
  // Scale by 1e8 (the precision the rest of ENDF uses for similar fields) and
  // decode digit-by-digit. llround absorbs the FP noise of values like 1.55
  // (1.5499999999999998 in IEEE 754).
  constexpr std::int64_t kScale = 100000000;  // 10^8
  const std::int64_t scaled = std::llround(rtyp * static_cast<double>(kScale));
  const std::int64_t ipart = scaled / kScale;
  seq.push(static_cast<int>(ipart));
  std::int64_t frac = scaled - ipart * kScale;
  if (frac == 0)
    return seq;
  // Count trailing zeros so 1.5 yields {1, 5} (not {1, 5, 0, 0, ...}) while
  // 1.04 still yields {1, 0, 4}.
  int trailingZeros = 0;
  for (std::int64_t f = frac; f % 10 == 0; f /= 10)
    ++trailingZeros;
  const int digitsToEmit = 8 - trailingZeros;
  std::int64_t divisor = kScale / 10;  // 10^7
  for (int k = 0; k < digitsToEmit; ++k) {
    const int d = static_cast<int>(frac / divisor);
    seq.push(d);
    frac -= static_cast<std::int64_t>(d) * divisor;
    divisor /= 10;
  }
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
