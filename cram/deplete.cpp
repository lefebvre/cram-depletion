#include "cram/deplete.hpp"

#include <stdexcept>
#include <utility>

namespace cram {
namespace {

constexpr double kBarn = 1e-24;  // cm^2 per barn
constexpr double kEvToJoule = 1.602176634e-19;

}  // namespace

ReactionXS reactionXs(const ChainReaction& channel, double sigmaBarn) {
  return ReactionXS{.type = channel.type,
                    .target = channel.target,
                    .sigma = sigmaBarn * channel.branching,
                    .q = channel.q};
}

DepletionSystem::DepletionSystem(const DepletionChain& chain)
    : chain_(chain), decay_(chain.decayMatrix()) {}

void DepletionSystem::requireFissionQ(const std::vector<ReactionXS>& reactions) {
  for (const auto& r : reactions) {
    if (r.type == ReactionType::Fission && !(r.q > 0.0))
      throw std::invalid_argument(
          "cram: a fission channel needs q > 0 under constant-power normalization");
  }
}

void DepletionSystem::setReactions(const Zai& parent, std::vector<ReactionXS> reactions) {
  for (const auto& r : reactions) {
    if (r.sigma < 0.0)
      throw std::invalid_argument("cram: negative cross section for " + parent.str());
  }
  if (norm_ == Normalization::ConstantPower && power_ != 0.0)
    requireFissionQ(reactions);
  for (auto& [p, rxns] : reactions_) {
    if (p == parent) {
      rxns = std::move(reactions);
      return;
    }
  }
  reactions_.emplace_back(parent, std::move(reactions));
}

void DepletionSystem::setConstantPower(double power) {
  if (power != 0.0) {
    for (const auto& [parent, rxns] : reactions_)
      requireFissionQ(rxns);
  }
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
  if (power_ == 0.0)
    return 0.0;
  const double w = fissionPowerWeight(n);
  if (!(w > 0.0))
    throw std::domain_error(
        "cram: constant-power normalization with no fissile material present "
        "(zero fission-power weight)");
  return power_ / w;
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
        } else if (r.target && chain_.indexOf(*r.target) >= 0) {
          chain_.addReaction(t, parent, *r.target, rate);
        } else {
          // Product not tracked in the chain (an untracked target, or one the
          // chain never registered): the parent is still consumed, the product
          // is just lost.
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
