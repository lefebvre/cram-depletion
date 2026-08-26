#pragma once
//
// A constant-cross-section depletion problem: a DepletionChain (decay data +
// fission yields + reaction topology) together with one-group microscopic
// cross sections and a flux/power normalization. This is the transport-free
// analogue of OpenMC's IndependentOperator: the neutron spectrum is frozen
// (the one-group cross sections are fixed), but the flux MAGNITUDE may float
// to hold a target fission power as the fissile inventory depletes.
//
// That power normalization is what makes the burnup matrix A composition-
// dependent (A scales with the flux, and the flux depends on the current
// inventory), so the higher-order integrators in integrator.hpp are exercised
// rather than collapsing to a single exact exp(A*dt).
//
// One system describes one material: one chain, one cross-section set, one
// composition vector, one flux. A model with several burnable regions holds one
// DepletionSystem (and one Integrator) per region and shares power between
// them itself; that bookkeeping is deliberately not done here.
//
#include <Eigen/SparseCore>
#include <optional>
#include <utility>
#include <vector>

#include "cram/chain.hpp"
#include "cram/integrator.hpp"
#include "cram/reaction.hpp"

namespace cram {

// One-group microscopic cross section for a single reaction channel.
//
// `type`, `target` and `sigma` have no default member initializers, for the
// reason DecayMode's members have none: a forgotten sigma is a channel that
// silently never fires, and a forgotten target is a channel that silently
// destroys atoms. -Wmissing-field-initializers therefore requires all three at
// every braced initialization. `target` is optional in *value*, not in
// presence: std::nullopt is the explicit statement that the product is not
// tracked (an OpenMC "Nothing" target), in which case the parent is consumed
// and nothing is produced. Fission ignores the target; its products come from
// the chain's yield tables.
//
// `q` and `energy` keep initializers because zero and thermal are their correct
// values when not stated: q is read only for fission (and a fission channel
// under constant-power normalization must supply it -- see setConstantPower),
// and 0.0253 eV is the conventional incident energy for the yield lookup.
struct ReactionXS {
  ReactionType type;
  std::optional<Zai> target;  // product; nullopt = consumed, product untracked
  double sigma;               // microscopic cross section [barn]
  double q = 0.0;             // fission energy release [eV] (fission only)
  double energy = 0.0253;     // incident neutron energy [eV] for the yield lookup
};

// The cross-section entry for a topology channel with the given microscopic
// cross section [barn] for the whole channel. The channel's branching scales
// sigma, and its Q value becomes ReactionXS::q, which is how a fission channel
// read from an OpenMC chain reaches constant-power normalization without the
// caller copying fields by hand.
ReactionXS reactionXs(const ChainReaction& channel, double sigmaBarn);

class DepletionSystem {
public:
  enum class Normalization { ConstantPower, ConstantFlux };

  // The chain must outlive the system, which holds a reference to it; its decay
  // matrix is cached at construction (decay is composition-independent). The
  // rvalue overload is deleted so a temporary chain cannot be bound: the
  // reference would dangle as soon as the constructor returned.
  explicit DepletionSystem(const DepletionChain& chain);
  explicit DepletionSystem(DepletionChain&&) = delete;

  // Set the one-group reaction cross sections for `parent` (replaces any prior
  // set). Parents and products are expected to already exist in the chain.
  // Throws std::invalid_argument for a negative sigma, and for a fission
  // channel with q <= 0 while a non-zero constant power is active.
  void setReactions(const Zai& parent, std::vector<ReactionXS> reactions);

  // Hold a target fission power: the flux scales each assemble() so that
  // sum_j n_j * sigma_f,j * phi * Q_j == power (Q in eV, converted to J).
  // `power` is a power density in the same volume basis as the number densities
  // (e.g. W/cm^3 when n is atom/cm^3). Throws std::invalid_argument if any
  // fission channel already set has q <= 0: such a channel would contribute
  // fissions but no power, so the normalization would be silently wrong.
  void setConstantPower(double power);

  // Hold a fixed scalar flux [n/cm^2/s]; A is then composition-independent.
  void setConstantFlux(double flux);

  Normalization normalization() const { return norm_; }

  // Scalar flux [n/cm^2/s] that assemble() would use for composition n. Under
  // constant power with a non-zero target, a composition with no fissile
  // material (zero fission-power weight) cannot produce that power, and this
  // throws std::domain_error rather than returning a flux of zero -- the
  // latter would let a mis-keyed fission cross section make an entire
  // irradiation vanish without a diagnostic. A zero power target is
  // legitimately zero flux.
  double fluxFor(const Eigen::VectorXd& n) const;

  // Total fission power for composition n at the flux fluxFor(n) (diagnostic).
  double powerFor(const Eigen::VectorXd& n) const;

  // Assemble the burnup matrix A(n) [1/s] = (cached decay) + (reactions at the
  // normalized flux).
  Eigen::SparseMatrix<double> assemble(const Eigen::VectorXd& n) const;

  // assemble() bound as a MatrixBuilder for the integrators. Captures `this`;
  // the system must outlive the integrator that holds the builder.
  MatrixBuilder matrixBuilder() const;

  const DepletionChain& chain() const { return chain_; }
  const std::vector<std::pair<Zai, std::vector<ReactionXS>>>& reactions() const {
    return reactions_;
  }

private:
  // Per-atom fission-power weight sum: sum_j n_j * sigma_f,j[barn]*1e-24 * Q_jJ.
  // Multiplying by flux gives fission power.
  double fissionPowerWeight(const Eigen::VectorXd& n) const;

  static void requireFissionQ(const std::vector<ReactionXS>& reactions);

  const DepletionChain& chain_;
  Eigen::SparseMatrix<double> decay_;  // cached, composition-independent
  std::vector<std::pair<Zai, std::vector<ReactionXS>>> reactions_;
  Normalization norm_ = Normalization::ConstantFlux;
  double power_ = 0.0;
  double flux_ = 0.0;
};

}  // namespace cram
