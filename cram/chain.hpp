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
//
// Deliberately without default member initializers: every member is a required
// input and zero is not an "empty" branch but a plausible-looking inert one
// (rtyp 0 decodes to a gamma, i.e. no Z/A change, and branching 0 contributes
// no rate). Leaving them out makes -Wmissing-field-initializers demand all four
// at each braced initialization, so a forgotten field is a build failure rather
// than a nuclide that silently never transmutes. Construct with braces only;
// never declare a bare `DecayMode m;`.
struct DecayMode {
  double rtyp;       // ENDF decay-mode code (see nuclide.hpp)
  double branching;  // branching ratio for this mode (sum over modes ~ 1)
  int finalState;    // RFS: isomeric state of the daughter
  bool isFission;    // spontaneous fission -> products from SFY table
};

// Per-nuclide decay information.
struct DecayData {
  double halfLife = 0.0;             // seconds; 0 or inf => stable
  double decayConstant = 0.0;        // ln(2)/halfLife; 0 if stable
  double gammaEnergyPerDecay = 0.0;  // average EM (gamma+X-ray) energy per decay [eV]
  // The explicit initializer is what keeps designated-initializer construction
  // (DecayData{.halfLife = ...}) free of -Wmissing-field-initializers: every
  // other member has one, so GCC would flag this member alone.
  std::vector<DecayMode> modes = {};  // empty if stable
};

// Independent fission yields for one parent at one incident energy.
// (Use ENDF MT454 = independent yields, NOT MT459 = cumulative yields, when
//  you model the full decay chain explicitly, otherwise you double count.)
struct FissionYields {
  double energy = 0.0;  // incident neutron energy [eV]; 0 => spontaneous
  std::vector<std::pair<Zai, double>> products = {};  // (product, yield per fission)
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

  // Register every reachable decay daughter not already in the chain, so matrix
  // assembly can never silently drop a daughter's production. Added daughters
  // carry no decay data (stable terminators). Returns the number added; 0 means
  // the chain was already closed. Idempotent.
  int close();

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
  // A yield set together with its products already resolved to matrix indices.
  //
  // addFissionYields() registers every product anyway, and add() hands back the
  // index while doing so, so the resolution is free at insert time. Doing it
  // there instead of during assembly removes a hash lookup per product from the
  // hottest loop in matrix building: a full NFY table runs to ~1000 products per
  // fissioning parent, across tens of parents, rebuilt for every depletion
  // region. Indices are stable because nuclides_ only ever grows and add()
  // never reorders it.
  struct YieldEntry {
    FissionYields yields;                          // as supplied; nearestYields() hands this back
    std::vector<std::pair<int, double>> resolved;  // (matrix index, yield per fission)
  };

  std::vector<Zai> nuclides_;
  std::unordered_map<std::int64_t, int> index_;
  std::unordered_map<std::int64_t, DecayData> decay_;
  std::unordered_map<std::int64_t, std::vector<YieldEntry>> nfy_;

  void decayTriplets(std::vector<Eigen::Triplet<double>>& t) const;

  const YieldEntry* nearestEntry(const Zai& parent, double energy) const;

  // Emit production triplets for the given fission-yield set, each weighted by
  // `rate` * yield. Caller is responsible for the parent-removal term.
  void emitFissionProducts(std::vector<Eigen::Triplet<double>>& t, int parentIndex,
                           const YieldEntry& entry, double rate) const;
};

}  // namespace cram
