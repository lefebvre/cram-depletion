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
#include <utility>
#include <vector>

#include "cram/nuclide.hpp"

namespace cram {

constexpr double kLn2 = std::numbers::ln2;

// One decay branch of a parent nuclide.
//
// The first four members deliberately have no default member initializers:
// every one is a required input and zero is not an "empty" branch but a
// plausible-looking inert one (rtyp 0 decodes to a gamma, i.e. no Z/A change,
// and branching 0 contributes no rate). Leaving them out makes
// -Wmissing-field-initializers demand all four at each braced initialization,
// so a forgotten field is a build failure rather than a nuclide that silently
// never transmutes. Construct with braces only; never declare a bare
// `DecayMode m;`.
//
// `daughter` keeps an initializer because "absent" is its correct default, not a
// placeholder: an unset daughter means "derive it from rtyp and finalState with
// applyDecay()", which is what ENDF decay data implies. It is set only by
// sources that name the product explicitly (an OpenMC depletion chain gives a
// target per decay mode) so that the matrix honors the source's topology even
// where it disagrees with the RTYP transition rules. Ignored when the mode is
// fission.
struct DecayMode {
  double rtyp;                                 // ENDF decay-mode code (see nuclide.hpp)
  double branching;                            // branching ratio for this mode (sum over modes ~ 1)
  int finalState;                              // RFS: isomeric state of the daughter
  bool isFission;                              // spontaneous fission -> products from SFY table
  std::optional<Zai> daughter = std::nullopt;  // explicit product; else derived
};

// Per-nuclide decay information.
//
// halfLife has no default member initializer, for the reason DecayMode's members
// have none: a half-life of 0 reads as "stable" everywhere downstream, so an
// omitted one produces a nuclide that is silently frozen in place rather than an
// obviously invalid one. -Wmissing-field-initializers therefore requires it at
// every braced initialization.
//
// The remaining members keep theirs because zero is the correct value, not a
// placeholder — and a member with an initializer is exempt from that warning, so
// each one below is a deliberate statement that it need not be supplied:
//   * decayConstant is derived, not input: setDecay() recomputes it from
//     halfLife and overwrites whatever the caller passed.
//   * gammaEnergyPerDecay is optional; MT457 omits it for many nuclides, and
//     0 eV of emitted EM energy is what "not reported" means here.
//   * modes is empty for a stable nuclide, and is otherwise appended to after
//     construction, so it cannot be required at the brace.
struct DecayData {
  double halfLife;                    // seconds; 0 or inf => stable
  double decayConstant = 0.0;         // ln(2)/halfLife; 0 if stable. Set by setDecay()
  double gammaEnergyPerDecay = 0.0;   // average EM (gamma+X-ray) energy per decay [eV]
  std::vector<DecayMode> modes = {};  // empty if stable
};

// Independent fission yields for one parent at one incident energy.
// (Use ENDF MT454 = independent yields, NOT MT459 = cumulative yields, when
//  you model the full decay chain explicitly, otherwise you double count.)
//
// `energy` has no default member initializer: 0 already means something —
// spontaneous fission — so an omitted incident energy does not fail, it
// reclassifies a neutron-induced yield table as an SF one and sends
// nearestYields() to the wrong table for every lookup. `products` keeps its
// initializer because the product list is filled in after construction.
struct FissionYields {
  double energy;                                      // incident neutron energy [eV]
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
  //
  // If `droppedDaughters` is non-null it receives the number of decay daughters
  // that were not registered in the chain, and whose production is therefore
  // absent from the matrix — meaning the matrix does not conserve atoms. Call
  // close() beforehand and it will be zero. When the pointer is null the same
  // condition is reported as a warning on stderr instead, so it can never pass
  // unnoticed either way.
  Eigen::SparseMatrix<double> decayMatrix(int* droppedDaughters = nullptr) const;

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

  // The single tracked product of a decay mode: the explicit `daughter` when
  // one was supplied, otherwise the nuclide applyDecay() derives from the RTYP
  // sequence. Empty when the mode is (or encodes) fission, whose products come
  // from the yield tables instead. Shared by close() and decayTriplets() so the
  // daughters registered by the one are exactly the daughters produced by the
  // other, and public so a caller contracting over the decay topology (the
  // sensitivity code does) sees the same product the matrix does.
  static std::optional<Zai> decayDaughter(const Zai& parent, const DecayMode& m);

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
  // Matrix indices of the nuclides carrying decay data, in the order setDecay()
  // first registered them. Matrix assembly walks this instead of decay_ so the
  // triplet order — and therefore the summation order of duplicates, and the
  // last bits of the result — cannot depend on hash iteration order. Kept in
  // sync in setDecay(), the only place decay_ gains a key.
  std::vector<int> decayOrder_;
  std::unordered_map<std::int64_t, std::vector<YieldEntry>> nfy_;

  // `dropped` receives the number of decay daughters missing from the chain.
  void decayTriplets(std::vector<Eigen::Triplet<double>>& t, int& dropped) const;

  static void warnDroppedDaughters(int dropped);

  const YieldEntry* nearestEntry(const Zai& parent, double energy) const;

  // Emit production triplets for the given fission-yield set, each weighted by
  // `rate` * yield. Caller is responsible for the parent-removal term.
  void emitFissionProducts(std::vector<Eigen::Triplet<double>>& t, int parentIndex,
                           const YieldEntry& entry, double rate) const;
};

}  // namespace cram
