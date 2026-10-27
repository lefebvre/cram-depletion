#include "cram/chain.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cram {

int DepletionChain::add(const Zai& z) {
  if (!z.representable())
    throw std::invalid_argument("cram: nuclide " + z.str() +
                                " is outside the range a ZAI key can hold "
                                "(Z >= 0, 1 <= A <= 999, 0 <= I <= 9)");
  auto it = index_.find(z.key());
  if (it != index_.end())
    return it->second;
  int idx = static_cast<int>(nuclides_.size());
  nuclides_.push_back(z);
  index_.emplace(z.key(), idx);
  ++revision_;
  return idx;
}

int DepletionChain::indexOf(const Zai& z) const {
  if (!z.representable())
    return -1;  // its key would alias a registered nuclide's
  auto it = index_.find(z.key());
  return it == index_.end() ? -1 : it->second;
}

void DepletionChain::setDecay(const Zai& z, DecayData d) {
  if (d.halfLife > 0.0 && std::isfinite(d.halfLife))
    d.decayConstant = kLn2 / d.halfLife;
  else
    d.decayConstant = 0.0;
  const int idx = add(z);
  // Default-construct on first insert, then assign unconditionally, so `d` is
  // moved from exactly once on every path. Passing std::move(d) to try_emplace
  // and moving again in the not-inserted branch would also be correct --
  // try_emplace is specified not to move from its arguments when no insertion
  // happens -- but it reads as a use-after-move and static analysis flags it as
  // one, which is not worth an inline suppression here.
  const auto [it, inserted] = decay_.try_emplace(z.key());
  if (inserted)
    decayOrder_.push_back(idx);  // see decayTriplets(): fixes assembly order
  it->second = std::move(d);
  ++revision_;
}

void DepletionChain::addFissionYields(const Zai& parent, FissionYields y) {
  add(parent);
  std::vector<std::pair<int, double>> resolved;
  resolved.reserve(y.products.size());
  // add() registers the product and returns its index in one step, so the
  // matrix index is known here and never has to be looked up during assembly.
  for (const auto& [p, yield] : y.products)
    resolved.emplace_back(add(p), yield);
  // Built in one brace rather than default-constructed and filled: YieldEntry
  // holds a FissionYields, whose `energy` has no default initializer, so a bare
  // `YieldEntry entry;` would leave it indeterminate until the assignment.
  YieldEntry entry{.yields = std::move(y), .resolved = std::move(resolved)};
  nfy_[parent.key()].push_back(std::move(entry));
  ++revision_;
}

const DecayData* DepletionChain::decay(const Zai& z) const {
  auto it = decay_.find(z.key());
  return it == decay_.end() ? nullptr : &it->second;
}

int DepletionChain::close() {
  // Snapshot the decaying-nuclide keys first: add() below mutates nuclides_/index_
  // (not decay_), and daughters added here carry no decay data, so one pass closes
  // the chain. Mirrors the daughter derivation in decayTriplets().
  std::vector<std::int64_t> parents;
  parents.reserve(decay_.size());
  for (const auto& [key, d] : decay_) {
    if (d.decayConstant != 0.0)
      parents.push_back(key);
  }
  // Sort before walking. decay_ is an unordered_map, so the order it yields
  // depends on the hash function and the rehash history -- which differ between
  // standard library implementations. Since add() assigns matrix indices in call
  // order, an unsorted walk would make the nuclide-to-index mapping, and with it
  // the matrix layout, the LU pivoting and the low-order bits of every result,
  // vary from platform to platform. Sorting costs O(n log n) once per close()
  // and buys a chain ordering that is a function of the data alone.
  std::sort(parents.begin(), parents.end());

  int added = 0;
  for (std::int64_t key : parents) {
    const Zai parent{.z = static_cast<int>(key / 10000),
                     .a = static_cast<int>((key / 10) % 1000),
                     .i = static_cast<int>(key % 10)};
    for (const DecayMode& m : decay_.at(key).modes) {
      const std::optional<Zai> daughter = decayDaughter(parent, m);
      if (!daughter || indexOf(*daughter) >= 0)
        continue;      // fission products come from the SFY/NFY tables, not one daughter
      add(*daughter);  // bare stable terminator -> production is never dropped
      ++added;
    }
  }
  return added;
}

const DepletionChain::YieldEntry* DepletionChain::nearestEntry(const Zai& parent,
                                                               double energy) const {
  // Energy 0 is the spontaneous-fission slot, asked for by name: it gets the
  // SF table or nothing, never the neutron-induced table nearest to 0 eV.
  if (energy == 0.0)
    return spontaneousFissionEntry(parent);
  auto it = nfy_.find(parent.key());
  if (it == nfy_.end() || it->second.empty())
    return nullptr;
  const YieldEntry* best = nullptr;
  double bestDist = std::numeric_limits<double>::infinity();
  for (const auto& entry : it->second) {
    if (entry.yields.energy == 0.0)
      continue;  // the spontaneous-fission table: not a neutron-induced one
    double d = std::abs(entry.yields.energy - energy);
    if (d < bestDist) {
      bestDist = d;
      best = &entry;
    }
  }
  return best;
}

const FissionYields* DepletionChain::nearestYields(const Zai& parent, double energy) const {
  const YieldEntry* entry = nearestEntry(parent, energy);
  return entry == nullptr ? nullptr : &entry->yields;
}

const DepletionChain::YieldEntry* DepletionChain::spontaneousFissionEntry(const Zai& parent) const {
  auto it = nfy_.find(parent.key());
  if (it == nfy_.end())
    return nullptr;
  const auto sf = std::find_if(it->second.begin(), it->second.end(),
                               [](const YieldEntry& e) { return e.yields.energy == 0.0; });
  return sf == it->second.end() ? nullptr : &*sf;
}

const FissionYields* DepletionChain::spontaneousFissionYields(const Zai& parent) const {
  const YieldEntry* entry = spontaneousFissionEntry(parent);
  return entry == nullptr ? nullptr : &entry->yields;
}

void DepletionChain::addReaction(std::vector<Eigen::Triplet<double>>& t, const Zai& parent,
                                 const Zai& product, double rate) const {
  int ip = indexOf(parent);
  int jp = indexOf(product);
  if (ip < 0 || jp < 0 || rate == 0.0)
    return;
  t.emplace_back(jp, ip, rate);   // production of product
  t.emplace_back(ip, ip, -rate);  // removal of parent
}

void DepletionChain::addFissionSource(std::vector<Eigen::Triplet<double>>& t, const Zai& parent,
                                      double fissionRate, double energy) const {
  int ip = indexOf(parent);
  if (ip < 0 || fissionRate == 0.0)
    return;
  const YieldEntry* entry = nearestEntry(parent, energy);
  if (entry == nullptr)
    return;
  t.emplace_back(ip, ip, -fissionRate);  // parent consumed by fission
  emitFissionProducts(t, ip, *entry, fissionRate);
}

void DepletionChain::emitFissionProducts(std::vector<Eigen::Triplet<double>>& t, int parentIndex,
                                         const YieldEntry& entry, double rate) const {
  // No index lookup and no validity check: addFissionYields() registered every
  // product and recorded its index, so all of them are in the chain.
  for (const auto& [j, yield] : entry.resolved)
    t.emplace_back(j, parentIndex, rate * yield);
}

void DepletionChain::decayTriplets(std::vector<Eigen::Triplet<double>>& t, int& dropped,
                                   int& untrackedSf) const {
  dropped = 0;      // daughters not registered -> production lost (chain not closed)
  untrackedSf = 0;  // SF modes with no SF table -> products lost
  // Walk decayOrder_, not decay_. Triplet order reaches setFromTriplets(), which
  // sums duplicate entries in the order it receives them, so iterating an
  // unordered_map would let the hash function decide the floating-point result.
  // decayOrder_ records the matrix index of each nuclide as setDecay() first
  // registered it, so this walk is a function of the input alone.
  //
  // It holds indices rather than being a scan over nuclides_ for cost reasons:
  // once close() has run, nuclides_ contains every stable daughter too -- more
  // than twice the entries in a realistic chain -- and scanning it measured
  // 15-21% slower than the hash-ordered walk it replaced. Indexing straight to
  // the decaying nuclides keeps the iteration count where it was.
  //
  // This still costs 5-14% against that hash-ordered walk, which had each
  // DecayData in hand as it iterated while this hops between decayOrder_,
  // nuclides_ and decay_. Caching DecayData pointers here would remove the
  // lookup, but assembly is well under 1% of a region's total cost, so the
  // residual is ~0.1% overall -- not worth pinning the code to the reference
  // stability of an unordered_map. Determinism is the point of this walk.
  std::vector<std::pair<int, double>> column;  // reused across nuclides
  for (const int i : decayOrder_) {
    const DecayData& d = decay_.find(nuclides_[static_cast<std::size_t>(i)].key())->second;
    if (d.decayConstant == 0.0)
      continue;
    decayStencil(i, d, column, dropped, untrackedSf);
    for (const auto& [j, c] : column)
      t.emplace_back(j, i, d.decayConstant * c);
  }
}

void DepletionChain::decayStencil(int i, const DecayData& d,
                                  std::vector<std::pair<int, double>>& out, int& dropped,
                                  int& untrackedSf) const {
  out.clear();
  if (d.decayConstant == 0.0)
    return;
  const Zai& parent = nuclides_[static_cast<std::size_t>(i)];
  out.emplace_back(i, -1.0);  // total removal of parent

  for (const auto& m : d.modes) {
    if (m.branching == 0.0)
      continue;

    if (m.isFission) {
      // Only the parent's own SF table. Substituting the nearest neutron-
      // induced one would put a thermal fission's products where a
      // spontaneous fission's belong, with nothing to say it happened.
      if (const YieldEntry* entry = spontaneousFissionEntry(parent)) {
        // No index lookup and no validity check: addFissionYields() registered
        // every product and recorded its index.
        for (const auto& [j, yield] : entry->resolved)
          out.emplace_back(j, m.branching * yield);
      } else {
        ++untrackedSf;  // removed with the rest of the decay, nothing produced
      }
      continue;
    }

    const std::optional<Zai> daughter = decayDaughter(parent, m);
    if (!daughter)
      continue;  // RTYP encodes fission (e.g. 1.6): no single product
    const int j = indexOf(*daughter);
    if (j >= 0)
      out.emplace_back(j, m.branching);
    else
      ++dropped;  // unregistered daughter; surfaced by decayMatrix() instead of silently lost
  }
}

void DepletionChain::decayStencil(int index, std::vector<std::pair<int, double>>& out) const {
  if (index < 0 || index >= size())
    throw std::out_of_range("cram: DepletionChain::decayStencil index out of range");
  out.clear();
  const DecayData* d = decay(nuclides_[static_cast<std::size_t>(index)]);
  if (d == nullptr)
    return;
  int dropped = 0;
  int untrackedSf = 0;
  decayStencil(index, *d, out, dropped, untrackedSf);
}

void DepletionChain::reactionStencil(const Zai& parent, bool isFission,
                                     const std::optional<Zai>& target, double energy,
                                     std::vector<std::pair<int, double>>& out) const {
  out.clear();
  const int p = indexOf(parent);
  if (p < 0)
    return;
  out.emplace_back(p, -1.0);  // the parent is consumed whatever becomes of the product
  if (isFission) {
    if (const YieldEntry* entry = nearestEntry(parent, energy))
      for (const auto& [j, yield] : entry->resolved)
        out.emplace_back(j, yield);
  } else if (target) {
    const int j = indexOf(*target);
    if (j >= 0)
      out.emplace_back(j, 1.0);
  }
}

std::optional<Zai> DepletionChain::decayDaughter(const Zai& parent, const DecayMode& m) {
  if (m.isFission)
    return std::nullopt;
  if (m.daughter)
    return m.daughter;
  bool fission = false;
  const Zai derived = applyDecay(parent, m.rtyp, m.finalState, fission);
  if (fission)
    return std::nullopt;
  return derived;
}

void DepletionChain::warnDroppedDaughters(int dropped) {
  // Do not fail silently: a non-zero count means the chain was not closed and the
  // decay matrix does not conserve atoms. Call DepletionChain::close() beforehand.
  if (dropped > 0) {
    std::fprintf(stderr,
                 "cram: WARNING - %d decay daughter(s) not registered; their production "
                 "is dropped and the decay matrix will not conserve atoms. "
                 "Call DepletionChain::close() before building the matrix.\n",
                 dropped);
  }
}

void DepletionChain::warnUntrackedSpontaneousFission(int untracked) {
  if (untracked > 0) {
    std::fprintf(stderr,
                 "cram: WARNING - %d spontaneous-fission mode(s) belong to a parent with no "
                 "spontaneous-fission yield table (energy 0); the parent is removed but no "
                 "products are produced, so the decay matrix will not conserve atoms.\n",
                 untracked);
  }
}

Eigen::SparseMatrix<double> DepletionChain::decayMatrix(int* droppedDaughters,
                                                        int* untrackedSpontaneousFission) const {
  std::vector<Eigen::Triplet<double>> t;
  // Each decaying nuclide contributes a diagonal term plus one production
  // triplet per decay mode (spontaneous-fission modes add more, but those are
  // a small minority — this is just a sizing hint).
  std::size_t estimate = 0;
  for (const auto& [_, d] : decay_) {
    if (d.decayConstant == 0.0)
      continue;
    estimate += 1 + d.modes.size();
  }
  t.reserve(estimate);

  int dropped = 0;
  int untrackedSf = 0;
  decayTriplets(t, dropped, untrackedSf);
  // Reported through the out-parameter when the caller asked for it, and only
  // written to stderr otherwise. A caller that inspects the count is handling
  // the condition and does not need a warning printed underneath it. Passing it
  // out rather than caching it on the object also keeps decayMatrix() free of
  // mutable state, so a chain shared read-only across region worker threads can
  // be assembled from concurrently.
  if (droppedDaughters != nullptr)
    *droppedDaughters = dropped;
  else
    warnDroppedDaughters(dropped);
  if (untrackedSpontaneousFission != nullptr)
    *untrackedSpontaneousFission = untrackedSf;
  else
    warnUntrackedSpontaneousFission(untrackedSf);

  return finalize(std::move(t));
}

Eigen::SparseMatrix<double> DepletionChain::finalize(
    std::vector<Eigen::Triplet<double>> triplets) const {
  Eigen::SparseMatrix<double> A(size(), size());
  A.setFromTriplets(triplets.begin(), triplets.end());  // sums duplicates
  A.makeCompressed();
  return A;
}

}  // namespace cram
