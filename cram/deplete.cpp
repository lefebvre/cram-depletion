#include "cram/deplete.hpp"

#include <utility>

namespace cram {
namespace {

constexpr double kBarn = 1e-24;  // cm^2 per barn
constexpr double kEvToJoule = 1.602176634e-19;

}  // namespace

Zai reactionProduct(const Zai& parent, ReactionType type) {
  Zai p{parent.z, parent.a, 0};
  switch (type) {
    case ReactionType::NGamma:  // (n,gamma): A -> A+1
      p.a += 1;
      break;
    case ReactionType::N2n:  // (n,2n): A -> A-1
      p.a -= 1;
      break;
    case ReactionType::N3n:  // (n,3n): A -> A-2
      p.a -= 2;
      break;
    case ReactionType::N4n:  // (n,4n): A -> A-3
      p.a -= 3;
      break;
    case ReactionType::NAlpha:  // (n,alpha): Z-2, A-3
      p.z -= 2;
      p.a -= 3;
      break;
    case ReactionType::NProton:  // (n,p): Z-1, A unchanged (A+1-1)
      p.z -= 1;
      break;
    case ReactionType::Fission:  // no single product
      break;
  }
  return p;
}

DepletionSystem::DepletionSystem(const DepletionChain& chain)
    : chain_(chain), decay_(chain.decayMatrix()) {}

void DepletionSystem::setReactions(const Zai& parent, std::vector<ReactionXS> reactions) {
  for (auto& [p, rxns] : reactions_) {
    if (p == parent) {
      rxns = std::move(reactions);
      return;
    }
  }
  reactions_.emplace_back(parent, std::move(reactions));
}

void DepletionSystem::setConstantPower(double power) {
  norm_ = Normalization::ConstantPower;
  power_ = power;
}

void DepletionSystem::setConstantFlux(double flux) {
  norm_ = Normalization::ConstantFlux;
  flux_ = flux;
}

double DepletionSystem::fissionPowerWeight(const Eigen::VectorXd& n) const {
  double w = 0.0;
  for (const auto& [parent, rxns] : reactions_) {
    int pi = chain_.indexOf(parent);
    if (pi < 0)
      continue;
    for (const auto& r : rxns) {
      if (r.type == ReactionType::Fission)
        w += n(pi) * r.sigma * kBarn * r.q * kEvToJoule;
    }
  }
  return w;
}

double DepletionSystem::fluxFor(const Eigen::VectorXd& n) const {
  if (norm_ == Normalization::ConstantFlux)
    return flux_;
  const double w = fissionPowerWeight(n);
  return w > 0.0 ? power_ / w : 0.0;
}

double DepletionSystem::powerFor(const Eigen::VectorXd& n) const {
  return fluxFor(n) * fissionPowerWeight(n);
}

Eigen::SparseMatrix<double> DepletionSystem::assemble(const Eigen::VectorXd& n) const {
  const double flux = fluxFor(n);

  std::vector<Eigen::Triplet<double>> t;
  if (flux != 0.0) {
    for (const auto& [parent, rxns] : reactions_) {
      const int pi = chain_.indexOf(parent);
      if (pi < 0)
        continue;
      for (const auto& r : rxns) {
        const double rate = r.sigma * kBarn * flux;  // per-atom rate [1/s]
        if (rate == 0.0)
          continue;
        if (r.type == ReactionType::Fission) {
          chain_.addFissionSource(t, parent, rate, r.energy);
        } else if (chain_.indexOf(r.target) >= 0) {
          chain_.addReaction(t, parent, r.target, rate);
        } else {
          // Product not tracked in the chain (e.g. an OpenMC "Nothing"
          // target): the parent is still consumed, the product is just lost.
          t.emplace_back(pi, pi, -rate);
        }
      }
    }
  }
  Eigen::SparseMatrix<double> reactions = chain_.finalize(std::move(t));
  return decay_ + reactions;
}

MatrixBuilder DepletionSystem::matrixBuilder() const {
  return [this](const Eigen::VectorXd& n) { return assemble(n); };
}

}  // namespace cram
