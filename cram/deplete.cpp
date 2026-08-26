#include "cram/deplete.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
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
    : chain_(chain), decay_(chain.decayMatrix()), chainRevision_(chain.revision()) {}

void DepletionSystem::refreshChain() {
  decay_ = chain_.decayMatrix();
  chainRevision_ = chain_.revision();
}

void DepletionSystem::requireFissionQ(const std::vector<ReactionXS>& reactions) {
  const bool missingQ = std::any_of(reactions.begin(), reactions.end(), [](const ReactionXS& r) {
    return r.type == ReactionType::Fission && !(r.q > 0.0);
  });
  if (missingQ)
    throw std::invalid_argument(
        "cram: a fission channel needs q > 0 under constant-power normalization");
}

void DepletionSystem::setReactions(const Zai& parent, std::vector<ReactionXS> reactions) {
  const bool negativeSigma = std::any_of(reactions.begin(), reactions.end(),
                                         [](const ReactionXS& r) { return r.sigma < 0.0; });
  if (negativeSigma)
    throw std::invalid_argument("cram: negative cross section for " + parent.str());
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
  if (!std::isfinite(power) || power < 0.0)
    throw std::invalid_argument("cram: constant power must be finite and non-negative");
  if (power != 0.0) {
    for (const auto& [parent, rxns] : reactions_)
      requireFissionQ(rxns);
  }
  norm_ = Normalization::ConstantPower;
  power_ = power;
}

void DepletionSystem::setConstantFlux(double flux) {
  if (!std::isfinite(flux) || flux < 0.0)
    throw std::invalid_argument("cram: constant flux must be finite and non-negative");
  norm_ = Normalization::ConstantFlux;
  flux_ = flux;
}

double DepletionSystem::fissionPowerWeight(const Eigen::VectorXd& n) const {
  double w = 0.0;
  for (const auto& [parent, rxns] : reactions_) {
    int pi = chain_.indexOf(parent);
    if (pi < 0)
      continue;
    // w is the init value, so the terms are summed in the same order a running
    // total would have visited them and the last bits do not move.
    w = std::accumulate(rxns.begin(), rxns.end(), w, [&](double acc, const ReactionXS& r) {
      if (r.type != ReactionType::Fission)
        return acc;
      return acc + n(pi) * r.sigma * kBarn * r.q * kEvToJoule;
    });
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
  // The decay half was cached against the chain as it stood at construction,
  // while the reaction half below is sized from the chain as it stands now. A
  // chain that grew in between -- close() registering the daughters the
  // constructor warned about is the ordinary way it happens -- leaves the two
  // describing different chains. Refuse rather than assert (debug build) or
  // silently truncate the sum to the smaller one (release build).
  if (chain_.revision() != chainRevision_)
    throw std::logic_error(
        "cram: the depletion chain changed after this DepletionSystem cached its decay "
        "matrix; call DepletionSystem::refreshChain() (or rebuild the system) after "
        "modifying the chain, e.g. after DepletionChain::close()");

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
          if (chain_.nearestYields(parent, r.energy) != nullptr) {
            chain_.addFissionSource(t, parent, rate, r.energy);
          } else {
            // Fissionable with no yield table (a chain trimmed of its yields,
            // or a nuclide whose yields were never supplied): the parent is
            // still consumed, exactly as for an untracked product below. Left
            // to addFissionSource() the whole channel would drop out -- no
            // removal either -- while fissionPowerWeight() still credits its
            // energy release, so constant-power normalization would silently
            // under-burn the material to hold a power the matrix does not
            // deliver.
            t.emplace_back(pi, pi, -rate);
          }
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
  Eigen::SparseMatrix<double> reactionMatrix = chain_.finalize(std::move(t));
  return decay_ + reactionMatrix;
}

MatrixBuilder DepletionSystem::matrixBuilder() const {
  return [this](const Eigen::VectorXd& n) { return assemble(n); };
}

}  // namespace cram
