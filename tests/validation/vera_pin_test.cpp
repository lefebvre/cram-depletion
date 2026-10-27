// Regression validation against OpenMC-generated VERA pin reference data.
//
// This replays a depletion that OpenMC ran with its first-order PREDICTOR
// integrator: at each step OpenMC froze the burnup matrix at the start-of-step
// reaction rates and applied exp(A*dt). Given the same per-step one-group cross
// sections, the same flux, and the same depletion chain, this engine's
// predictor march must reproduce OpenMC's number densities to ~CRAM precision.
// It therefore checks our matrix assembly + CRAM, isolated from the transport /
// cross-section error that we do not control.
//
// The reference data is produced offline by validation/openmc/generate_vera_pin.py
// (OpenMC + a multi-GB nuclear-data library, not a build/test dependency). When
// the data directory is absent the test SKIPS, so a clean checkout still passes.
//
// Data schema (CSV, one directory per case under CRAM_VALIDATION_DATA_DIR):
//   chain.xml      OpenMC depletion_chain (decay, fission yields, reactions)
//   schedule.csv   step,dt_seconds,flux               (step 1..N; flux n/cm^2/s)
//   micro_xs.csv   nuclide,reaction,xs_barn           (fixed, applied every step)
//   density.csv    step,nuclide,atoms                 (step 0..N; step 0 = initial)
//
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cram/chain_xml.hpp"
#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "cram/nuclide.hpp"
#include "cram/reaction.hpp"
#include "validation/vera_case.hpp"

#ifndef CRAM_VALIDATION_DATA_DIR
#define CRAM_VALIDATION_DATA_DIR ""
#endif

using namespace cram;
namespace fs = std::filesystem;

namespace {

using cram_validation::VeraCase;

}  // namespace

// See the comparison below for how these were chosen.
constexpr double kRelTol = 1e-10;
constexpr double kAbsTol = 1e-12;

class VeraPin : public ::testing::TestWithParam<std::string> {};

TEST_P(VeraPin, ReproducesOpenMCPredictorTrajectory) {
  const fs::path root = fs::path(CRAM_VALIDATION_DATA_DIR) / GetParam();
  if (CRAM_VALIDATION_DATA_DIR[0] == '\0' || !fs::exists(root / "chain.xml"))
    GTEST_SKIP() << "no OpenMC reference data at " << root
                 << " (run validation/openmc/generate_vera_pin.py)";

  const VeraCase c = cram_validation::loadVeraCase(root);
  EXPECT_EQ(c.diagnostics.unparsedNuclides, 0);
  EXPECT_EQ(c.diagnostics.unparsedDecayTargets, 0);
  EXPECT_EQ(c.diagnostics.unparsedYieldProducts, 0);
  // A reaction absent from the chain for a nuclide is ignored by the loader
  // (OpenMC only applies reactions the chain defines), even if the micro-xs file
  // carries a tally for it.
  ASSERT_GT(c.steps, 0);

  // Build the system once from the fixed cross sections, then march the
  // predictor step by step, re-normalizing only the flux each step.
  DepletionSystem sys(c.chain);
  int fissionChannels = 0;
  for (const auto& [parent, rxns] : c.xs) {
    for (const ReactionXS& r : rxns) {
      if (r.type != ReactionType::Fission)
        continue;
      ++fissionChannels;
      EXPECT_GT(r.q, 0.0) << parent.str() << ": fission Q did not reach the cross section";
    }
    sys.setReactions(parent, rxns);
  }
  EXPECT_GT(fissionChannels, 0);
  auto integ = makeIntegrator(IntegratorKind::Predictor, sys.matrixBuilder());

  std::vector<Eigen::VectorXd> traj{c.refDensity[0]};
  for (int k = 1; k <= c.steps; ++k) {
    const auto s = static_cast<std::size_t>(k - 1);
    sys.setConstantFlux(c.flux[s]);
    traj.push_back(integ->step(traj.back(), c.dt[s]));
  }

  // Compare EVERY nuclide the reference reports, at EVERY step -- not a
  // curated subset, and not only end of life.
  //
  // A curated list used to stand here, justified as leaving out deep-chain
  // trace actinides "sensitive to metastable/branching subtleties of the
  // simplified chain". What it actually left out were the nuclides exposing two
  // defects in the chain reader: fission yields delegated with
  // <neutron_fission_yields parent="..."/>, and the target OpenMC names on an
  // `sf` decay mode. Read correctly, all of them agree to CRAM precision.
  //
  // So the tolerance is set near that precision rather than at the 1e-3 it
  // used to be, which let ten orders of magnitude of regression through: a
  // nuclide passes when |cram - openmc| <= kRelTol * openmc + kAbsTol * scale,
  // scale being the step's largest inventory. The absolute term is CRAM's
  // floor -- an error of ~1e-15 of the largest entry, which a trace nuclide
  // many orders below it cannot resolve relatively -- and it is also what a
  // nuclide the reference reports as exactly zero is held to, so the engine
  // cannot produce one the reference does not.
  //
  // Measured on vera_pin1a over all 14 steps: the worst relative error of any
  // nuclide above 1e-9 of the largest inventory is 1.1e-14, and the worst
  // error of any nuclide is 3.0e-15 of the largest inventory. kRelTol and
  // kAbsTol leave four and two orders of margin for other toolchains and Eigen
  // versions, and still fail on anything that is not round-off.
  const auto name = [&](int i) {
    const Zai& z = c.chain.nuclides()[static_cast<std::size_t>(i)];
    return elementSymbol(z.z) + std::to_string(z.a) + (z.i != 0 ? "_m" + std::to_string(z.i) : "");
  };
  int compared = 0;
  double worstRel = 0.0;
  double worstAbs = 0.0;
  std::string worstNuclide = "none";
  for (int k = 1; k <= c.steps; ++k) {
    const auto ks = static_cast<std::size_t>(k);
    const Eigen::VectorXd& ref = c.refDensity[ks];
    const Eigen::VectorXd& n = traj[ks];
    const double scale = ref.maxCoeff();
    for (int i = 0; i < ref.size(); ++i) {
      if (c.reported[static_cast<std::size_t>(i)] == 0)
        continue;  // the reference says nothing about this nuclide
      const double diff = std::abs(n(i) - ref(i));
      worstAbs = std::max(worstAbs, diff / scale);
      if (ref(i) > 0.0) {
        ++compared;
        if (diff / ref(i) > worstRel && ref(i) > 1e-9 * scale) {  // above trace level
          worstRel = diff / ref(i);
          worstNuclide = name(i) + " step " + std::to_string(k);
        }
      }
      EXPECT_LE(diff, kRelTol * ref(i) + kAbsTol * scale)
          << name(i) << " step " << k << ": cram=" << n(i) << " openmc=" << ref(i);
    }
  }
  EXPECT_GT(compared, 50 * c.steps) << "the reference reports too few nuclides to be the pin case";
  RecordProperty("compared", compared);
  std::ostringstream text;  // std::to_string would print 4e-13 as 0.000000
  text << std::scientific << std::setprecision(3) << worstRel;
  RecordProperty("worst_relative_error", text.str());
  RecordProperty("worst_nuclide", worstNuclide);
  text.str("");
  text << worstAbs;
  RecordProperty("worst_error_relative_to_largest_inventory", text.str());

  const Eigen::VectorXd& ref = c.refDensity[static_cast<std::size_t>(c.steps)];
  // Every nuclide the benchmark reports must be among them: the comparison
  // above skips a nuclide the reference does not carry, so without this a data
  // file missing half the actinides would still pass it.
  for (const char* wanted : cram_validation::kBenchmarkNuclides) {
    const std::optional<Zai> z = parseNuclideName(wanted);
    ASSERT_TRUE(z.has_value()) << wanted;
    const int i = c.chain.indexOf(*z);
    EXPECT_TRUE(i >= 0 && ref(i) > 0.0) << wanted << ": absent from the reference data";
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, VeraPin, ::testing::Values("vera_pin1a"));
