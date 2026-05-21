#include "cram/chain.hpp"

#include <algorithm>
#include <cmath>

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
  if (nearestYields(parent, energy) == nullptr)
    return;
  t.emplace_back(ip, ip, -fissionRate);  // parent consumed by fission
  emitFissionProducts(t, ip, parent, fissionRate, energy);
}

void DepletionChain::emitFissionProducts(std::vector<Eigen::Triplet<double>>& t, int parentIndex,
                                         const Zai& parent, double rate, double energy) const {
  const FissionYields* y = nearestYields(parent, energy);
  if (y == nullptr)
    return;
  for (const auto& [prod, yield] : y->products) {
    int j = indexOf(prod);
    if (j >= 0)
      t.emplace_back(j, parentIndex, rate * yield);
  }
}

void DepletionChain::decayTriplets(std::vector<Eigen::Triplet<double>>& t) const {
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
        emitFissionProducts(t, i, parent, rate, 0.0);  // spontaneous fission
        continue;
      }

      bool fission = false;
      Zai daughter = applyDecay(parent, m.rtyp, m.finalState, fission);
      if (fission)
        continue;
      int j = indexOf(daughter);
      if (j >= 0)
        t.emplace_back(j, i, rate);
      // If j < 0 the daughter wasn't registered; production is dropped.
      // Register all reachable daughters when building the chain to avoid this.
    }
  }
}

Eigen::SparseMatrix<double> DepletionChain::decayMatrix() const {
  std::vector<Eigen::Triplet<double>> t;
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
