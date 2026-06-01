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
#include <Eigen/SparseCore>
#include <vector>

#include "cram/chain.hpp"
#include "cram/integrator.hpp"

namespace cram {

// Neutron reactions handled here. The product nuclide is stored explicitly on
// each ReactionXS (matching OpenMC's chain, where every reaction names its
// target); reactionProduct() computes the conventional target for a given
// parent when one is not supplied.
enum class ReactionType { Fission, NGamma, N2n, N3n, N4n, NAlpha, NProton };

// Conventional product of `type` on `parent` (ground state). Fission has no
// single product (it returns `parent` unchanged; products come from the NFY
// table) — callers should branch on the type instead.
Zai reactionProduct(const Zai& parent, ReactionType type);

// One-group microscopic cross section for a single reaction channel.
struct ReactionXS {
  ReactionType type = ReactionType::NGamma;
  Zai target;              // product (ignored for fission)
  double sigma = 0.0;      // microscopic cross section [barn]
  double q = 0.0;          // fission energy release [eV] (fission only)
  double energy = 0.0253;  // incident neutron energy [eV] for NFY lookup (fission)
};

class DepletionSystem {
public:
  enum class Normalization { ConstantPower, ConstantFlux };

  // The chain must outlive the system; its decay matrix is cached at
  // construction (decay is composition-independent).
  explicit DepletionSystem(const DepletionChain& chain);

  // Set the one-group reaction cross sections for `parent` (replaces any prior
  // set). Parents and products are expected to already exist in the chain.
  void setReactions(const Zai& parent, std::vector<ReactionXS> reactions);

  // Hold a target fission power: the flux scales each assemble() so that
  // sum_j n_j * sigma_f,j * phi * Q_j == power (Q in eV, converted to J). With
  // no fissile present the flux is zero. `power` is a power density in the same
  // volume basis as the number densities (e.g. W/cm^3 when n is atom/cm^3).
  void setConstantPower(double power);

  // Hold a fixed scalar flux [n/cm^2/s]; A is then composition-independent.
  void setConstantFlux(double flux);

  // Scalar flux [n/cm^2/s] that assemble() would use for composition n.
  double fluxFor(const Eigen::VectorXd& n) const;

  // Total fission power for composition n at the flux fluxFor(n) (diagnostic).
  double powerFor(const Eigen::VectorXd& n) const;

  // Assemble the burnup matrix A(n) [1/s] = (cached decay) + (reactions at the
  // normalized flux).
  Eigen::SparseMatrix<double> assemble(const Eigen::VectorXd& n) const;

  // assemble() bound as a MatrixBuilder for the integrators.
  MatrixBuilder matrixBuilder() const;

  const DepletionChain& chain() const { return chain_; }

private:
  // Per-atom fission-power weight sum: sum_j n_j * sigma_f,j[barn]*1e-24 * Q_jJ.
  // Multiplying by flux gives fission power.
  double fissionPowerWeight(const Eigen::VectorXd& n) const;

  const DepletionChain& chain_;
  Eigen::SparseMatrix<double> decay_;  // cached, composition-independent
  std::vector<std::pair<Zai, std::vector<ReactionXS>>> reactions_;
  Normalization norm_ = Normalization::ConstantFlux;
  double power_ = 0.0;
  double flux_ = 0.0;
};

}  // namespace cram
