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
  add(z);
  decay_[z.key()] = std::move(d);
}

void DepletionChain::addFissionYields(const Zai& parent, FissionYields y) {
  add(parent);
  for (auto& [p, _] : y.products)
    add(p);
  nfy_[parent.key()].push_back(std::move(y));
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
  int added = 0;
  for (std::int64_t key : parents) {
    Zai parent;
    parent.z = static_cast<int>(key / 10000);
    parent.a = static_cast<int>((key / 10) % 1000);
    parent.i = static_cast<int>(key % 10);
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

const FissionYields* DepletionChain::nearestYields(const Zai& parent, double energy) const {
  auto it = nfy_.find(parent.key());
  if (it == nfy_.end() || it->second.empty())
    return nullptr;
  const FissionYields* best = nullptr;
  double bestDist = std::numeric_limits<double>::infinity();
  for (const auto& y : it->second) {
    double d = std::abs(y.energy - energy);
    if (d < bestDist) {
      bestDist = d;
      best = &y;
    }
  }
  return best;
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
  const FissionYields* y = nearestYields(parent, energy);
  if (y == nullptr)
    return;
  t.emplace_back(ip, ip, -fissionRate);  // parent consumed by fission
  emitFissionProducts(t, ip, *y, fissionRate);
}

void DepletionChain::emitFissionProducts(std::vector<Eigen::Triplet<double>>& t, int parentIndex,
                                         const FissionYields& yields, double rate) const {
  for (const auto& [prod, yield] : yields.products) {
    int j = indexOf(prod);
    if (j >= 0)
      t.emplace_back(j, parentIndex, rate * yield);
  }
}

void DepletionChain::decayTriplets(std::vector<Eigen::Triplet<double>>& t) const {
  int dropped = 0;  // daughters not registered -> production lost (chain not closed)
  for (const auto& [key, d] : decay_) {
    if (d.decayConstant == 0.0)
      continue;
    int i = index_.at(key);
    Zai parent;
    parent.z = static_cast<int>(key / 10000);
    parent.a = static_cast<int>((key / 10) % 1000);
    parent.i = static_cast<int>(key % 10);

    t.emplace_back(i, i, -d.decayConstant);  // total removal of parent

    for (const auto& m : d.modes) {
      double rate = d.decayConstant * m.branching;
      if (rate == 0.0)
        continue;

      if (m.isFission) {
        if (const FissionYields* y = nearestYields(parent, 0.0))
          emitFissionProducts(t, i, *y, rate);  // spontaneous fission
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

Eigen::SparseMatrix<double> DepletionChain::decayMatrix() const {
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
  decayTriplets(t);
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
