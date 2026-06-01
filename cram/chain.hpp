#pragma once
//
// The depletion chain: the set of nuclides, their decay data, and fission
// yields, plus assembly of the burnup (transmutation) matrix A used by CRAM.
//
// Sign / layout convention (matches OpenMC and Pusa):
//   A is N x N, A(j, i) = rate at which nuclide i produces nuclide j.
//   The diagonal A(i, i) is negative and equals the total removal rate of i.
//   dn/dt = A n  ,  n(t) = exp(A t) n(0).
//
#include <Eigen/SparseCore>
#include <cstdint>
#include <numbers>
#include <optional>
#include <unordered_map>
#include <vector>

#include "cram/nuclide.hpp"

namespace cram {

constexpr double kLn2 = std::numbers::ln2;

// One decay branch of a parent nuclide.
struct DecayMode {
  double rtyp = 0.0;       // ENDF decay-mode code (see nuclide.hpp)
  double branching = 0.0;  // branching ratio for this mode (sum over modes ~ 1)
  int finalState = 0;      // RFS: isomeric state of the daughter
  bool isFission = false;  // spontaneous fission -> products from SFY table
  // Optional explicit daughter. When hasDaughter is true (and isFission is
  // false) the daughter below is used directly instead of deriving it from
  // rtyp via applyDecay(). This lets readers that already carry an explicit
  // target (e.g. an OpenMC depletion_chain XML) honor it exactly.
  bool hasDaughter = false;
  Zai daughter{};
};

// Per-nuclide decay information.
struct DecayData {
  double halfLife = 0.0;         // seconds; 0 or inf => stable
  double decayConstant = 0.0;    // ln(2)/halfLife; 0 if stable
  std::vector<DecayMode> modes;  // empty if stable
};

// Independent fission yields for one parent at one incident energy.
// (Use ENDF MT454 = independent yields, NOT MT459 = cumulative yields, when
//  you model the full decay chain explicitly, otherwise you double count.)
struct FissionYields {
  double energy = 0.0;                           // incident neutron energy [eV]; 0 => spontaneous
  std::vector<std::pair<Zai, double>> products;  // (product, yield per fission)
};

class DepletionChain {
public:
  // Register a nuclide (idempotent) and return its matrix index.
  int add(const Zai& z);

  // Look up an existing nuclide; returns -1 if absent.
  int indexOf(const Zai& z) const;

  int size() const { return static_cast<int>(nuclides_.size()); }
  const std::vector<Zai>& nuclides() const { return nuclides_; }

  void setDecay(const Zai& z, DecayData d);
  void addFissionYields(const Zai& parent, FissionYields y);

  const DecayData* decay(const Zai& z) const;
  // Independent yields for `parent` at the energy nearest to `energy` [eV].
  const FissionYields* nearestYields(const Zai& parent, double energy) const;

  // --- Matrix assembly -----------------------------------------------------

  // Pure radioactive decay (+ spontaneous fission if SFY were supplied as
  // yields at energy 0). This is all you need for a decay-only calculation.
  Eigen::SparseMatrix<double> decayMatrix() const;

  // Add a neutron-induced fission source for one parent into an existing
  // matrix builder. `fissionRate` is the per-atom fission rate of the parent
  // (= microscopic fission xs * scalar flux, units 1/s), evaluated at
  // `energy`. Removal of the parent and production of every fission product
  // (weighted by its independent yield) are added.
  void addFissionSource(std::vector<Eigen::Triplet<double>>& triplets, const Zai& parent,
                        double fissionRate, double energy) const;

  // Generic first-order removal/production term: A(prod,parent) += rate,
  // A(parent,parent) -= rate. Use for (n,gamma), (n,2n), etc.
  void addReaction(std::vector<Eigen::Triplet<double>>& triplets, const Zai& parent,
                   const Zai& product, double rate) const;

  // Build a sparse matrix from accumulated triplets (sums duplicates).
  Eigen::SparseMatrix<double> finalize(std::vector<Eigen::Triplet<double>> triplets) const;

private:
  std::vector<Zai> nuclides_;
  std::unordered_map<std::int64_t, int> index_;
  std::unordered_map<std::int64_t, DecayData> decay_;
  std::unordered_map<std::int64_t, std::vector<FissionYields>> nfy_;

  void decayTriplets(std::vector<Eigen::Triplet<double>>& t) const;

  // Emit production triplet(s) for one decay mode of `parent` (index `i`) at
  // the given branch rate. The parent-removal diagonal is the caller's job.
  void emitDecayMode(std::vector<Eigen::Triplet<double>>& t, const Zai& parent, int i,
                     const DecayMode& m, double rate) const;

  // Emit production triplets for the given fission-yield set, each weighted by
  // `rate` * yield. Caller is responsible for the parent-removal term.
  void emitFissionProducts(std::vector<Eigen::Triplet<double>>& t, int parentIndex,
                           const FissionYields& yields, double rate) const;
};

}  // namespace cram
