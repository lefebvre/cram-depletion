#include "cram/chain.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

namespace cram {

int DepletionChain::add(const Zai& z) {
  auto it = index_.find(z.key());
  if (it != index_.end())
    return it->second;
  int idx = static_cast<int>(nuclides_.size());
  nuclides_.push_back(z);
  index_.emplace(z.key(), idx);
  return idx;
}

int DepletionChain::indexOf(const Zai& z) const {
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
      if (m.isFission)
        continue;  // fission products come from the SFY/NFY tables, not one daughter
      bool fission = false;
      Zai daughter = applyDecay(parent, m.rtyp, m.finalState, fission);
      if (fission || indexOf(daughter) >= 0)
        continue;
      add(daughter);  // bare stable terminator -> production is never dropped
      ++added;
    }
  }
  return added;
}

const DepletionChain::YieldEntry* DepletionChain::nearestEntry(const Zai& parent,
                                                               double energy) const {
  auto it = nfy_.find(parent.key());
  if (it == nfy_.end() || it->second.empty())
    return nullptr;
  const YieldEntry* best = nullptr;
  double bestDist = std::numeric_limits<double>::infinity();
  for (const auto& entry : it->second) {
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

void DepletionChain::decayTriplets(std::vector<Eigen::Triplet<double>>& t, int& dropped) const {
  dropped = 0;  // daughters not registered -> production lost (chain not closed)
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
  for (const int i : decayOrder_) {
    const Zai& parent = nuclides_[static_cast<std::size_t>(i)];
    const DecayData& d = decay_.find(parent.key())->second;  // present by construction
    if (d.decayConstant == 0.0)
      continue;

    t.emplace_back(i, i, -d.decayConstant);  // total removal of parent

    for (const auto& m : d.modes) {
      double rate = d.decayConstant * m.branching;
      if (rate == 0.0)
        continue;

      if (m.isFission) {
        if (const YieldEntry* entry = nearestEntry(parent, 0.0))
          emitFissionProducts(t, i, *entry, rate);  // spontaneous fission
        continue;
      }

      bool fission = false;
      Zai daughter = applyDecay(parent, m.rtyp, m.finalState, fission);
      if (fission)
        continue;
      int j = indexOf(daughter);
      if (j >= 0)
        t.emplace_back(j, i, rate);
      else
        ++dropped;  // unregistered daughter; surfaced below instead of silently lost
    }
  }
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

Eigen::SparseMatrix<double> DepletionChain::decayMatrix(int* droppedDaughters) const {
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
  decayTriplets(t, dropped);
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
