#pragma once
//
// Loader for an OpenMC-generated depletion reference case (the data
// validation/openmc/generate_vera_pin.py writes). Shared by the regression test
// in tests/validation/vera_pin_test.cpp and the report generator in
// validation/report/, so the bridge from "OpenMC chain topology + a micro-xs
// table" to "a configured DepletionSystem" exists exactly once.
//
// Data schema (CSV, one directory per case):
//   chain.xml      OpenMC depletion_chain (decay, fission yields, reactions)
//   schedule.csv   step,dt_seconds,flux               (step 1..N; flux n/cm^2/s)
//   micro_xs.csv   nuclide,reaction,xs_barn           (fixed, applied every step)
//   density.csv    step,nuclide,atoms                 (step 0..N; step 0 = initial)
//
// Needs CRAM_WITH_CHAIN_XML: the chain arrives as OpenMC XML.
//
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "cram/chain_xml.hpp"
#include "cram/deplete.hpp"
#include "cram/integrator.hpp"
#include "cram/nuclide.hpp"
#include "cram/reaction.hpp"

namespace cram_validation {

namespace fs = std::filesystem;

// Split one CSV line, honoring double-quoted fields (reaction names like
// "(n,gamma)" contain commas). Quotes are stripped from the returned cells.
inline std::vector<std::string> splitCsv(const std::string& line) {
  std::vector<std::string> cols;
  std::string cell;
  bool inQuotes = false;
  for (char ch : line) {
    if (ch == '"') {
      inQuotes = !inQuotes;
    } else if (ch == ',' && !inQuotes) {
      cols.push_back(cell);
      cell.clear();
    } else {
      cell.push_back(ch);
    }
  }
  cols.push_back(cell);
  return cols;
}

// Rows of `p` with at least `columns` cells, header line skipped. A row with
// fewer is dropped rather than indexed past the end: these files are generated,
// but a truncated one is a plausible accident (an interrupted writer, a partial
// copy) and every caller below reads fixed column positions out of the result.
inline std::vector<std::vector<std::string>> readCsv(const fs::path& p, std::size_t columns) {
  std::vector<std::vector<std::string>> rows;
  std::ifstream in(p);
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    if (header) {  // skip the column-name row
      header = false;
      continue;
    }
    std::vector<std::string> cells = splitCsv(line);
    if (cells.size() >= columns)
      rows.push_back(std::move(cells));
  }
  return rows;
}

// One depletion case loaded from a data directory.
struct VeraCase {
  cram::DepletionChain chain;
  cram::ChainXmlDiagnostics diagnostics;
  std::vector<cram::ChainReaction> topology;
  // (parent, reaction) -> the chain's branch(es) for that reaction.
  std::map<std::pair<std::int64_t, cram::ReactionType>, std::vector<cram::ChainReaction>> branches;
  std::vector<double> dt;    // per step (s)
  std::vector<double> flux;  // per step (n/cm^2/s)
  std::unordered_map<cram::Zai, std::vector<cram::ReactionXS>, cram::ZaiHash> xs;
  std::vector<Eigen::VectorXd> refDensity;  // step 0..N
  int steps = 0;

  double totalSeconds() const {
    double t = 0.0;
    for (double d : dt)
      t += d;
    return t;
  }
};

// True when `root` holds a case this loader can read.
inline bool caseExists(const fs::path& root) {
  return fs::exists(root / "chain.xml") && fs::exists(root / "schedule.csv") &&
         fs::exists(root / "micro_xs.csv") && fs::exists(root / "density.csv");
}

// Load chain topology, schedule, cross sections and reference densities.
// Throws std::runtime_error (from the XML reader) if the chain cannot be read.
inline VeraCase loadVeraCase(const fs::path& root) {
  using namespace cram;
  VeraCase c;
  c.topology = loadDepletionChainXml(c.chain, (root / "chain.xml").string(), &c.diagnostics);
  for (const auto& r : c.topology)
    c.branches[{r.parent.key(), r.type}].push_back(r);

  // schedule.csv: step,dt_seconds,flux
  for (const auto& row : readCsv(root / "schedule.csv", 3)) {
    c.dt.push_back(std::stod(row[1]));
    c.flux.push_back(std::stod(row[2]));
  }
  c.steps = static_cast<int>(c.dt.size());

  // micro_xs.csv: nuclide,reaction,xs_barn. Only reactions the chain defines for
  // the nuclide are applied -- OpenMC ignores a tally the chain has no channel
  // for. reactionXs() splits a branched channel by branch fraction and carries
  // the fission Q through to constant-power normalization.
  for (const auto& row : readCsv(root / "micro_xs.csv", 3)) {
    const std::optional<Zai> parent = parseNuclideName(row[0]);
    const std::optional<ReactionType> type = reactionTypeFromName(row[1]);
    if (!parent || !type)
      continue;
    auto it = c.branches.find({parent->key(), *type});
    if (it == c.branches.end())
      continue;  // reaction not in the chain for this nuclide -> ignore the tally
    const double sigma = std::stod(row[2]);
    for (const ChainReaction& br : it->second)
      c.xs[*parent].push_back(reactionXs(br, sigma));
  }

  // density.csv: step,nuclide,atoms  (step 0..N)
  c.refDensity.assign(static_cast<std::size_t>(c.steps) + 1, Eigen::VectorXd::Zero(c.chain.size()));
  for (const auto& row : readCsv(root / "density.csv", 3)) {
    const int step = std::stoi(row[0]);
    const std::optional<Zai> z = parseNuclideName(row[1]);
    if (!z || step < 0 || step > c.steps)
      continue;
    const int idx = c.chain.indexOf(*z);
    if (idx >= 0)
      c.refDensity[static_cast<std::size_t>(step)](idx) = std::stod(row[2]);
  }
  return c;
}

// Configure a system from the case's fixed one-group cross sections.
inline void configure(const VeraCase& c, cram::DepletionSystem& sys) {
  for (const auto& [parent, rxns] : c.xs)
    sys.setReactions(parent, rxns);
}

// March the case with OpenMC's first-order predictor, re-normalizing only the
// flux each step. Returns the composition at step 0..N (step 0 is the reference
// initial condition), so a caller can compare the whole trajectory and not just
// end of life.
inline std::vector<Eigen::VectorXd> marchPredictor(const VeraCase& c) {
  cram::DepletionSystem sys(c.chain);
  configure(c, sys);
  auto integ = cram::makeIntegrator(cram::IntegratorKind::Predictor, sys.matrixBuilder());

  std::vector<Eigen::VectorXd> traj;
  traj.reserve(static_cast<std::size_t>(c.steps) + 1);
  Eigen::VectorXd n = c.refDensity[0];
  traj.push_back(n);
  for (int k = 1; k <= c.steps; ++k) {
    const auto s = static_cast<std::size_t>(k - 1);
    sys.setConstantFlux(c.flux[s]);
    n = integ->step(n, c.dt[s]);
    traj.push_back(n);
  }
  return traj;
}

// The benchmark's nuclides of interest: the actinides and fission products the
// VERA depletion benchmark reports (Yu & Forget, Ann. Nucl. Energy 170 (2022)
// 108973). Shared so the regression test and the validation report agree on
// what "a benchmark nuclide" means.
inline constexpr const char* kBenchmarkNuclides[] = {
    "U234",  "U235",  "U236",  "U238",  "Np237", "Pu238", "Pu239", "Pu240", "Pu241",
    "Pu242", "Am241", "Am243", "Xe135", "Cs137", "Nd148", "Sm149", "Gd157"};

// --- structural facts about the chain, for the report -----------------------
//
// These are properties of the input data, not judgements about the engine: they
// explain which nuclides the topology leaves under-determined.

// Nuclides carrying a fission channel for which the chain holds no usable yield
// table, so the fission products are unavailable. The parent is still consumed
// at the fission rate (see DepletionSystem::assemble); it is the products that
// the topology cannot place.
inline std::vector<cram::Zai> fissionWithoutYields(const VeraCase& c) {
  std::vector<cram::Zai> out;
  for (const auto& [key, brs] : c.branches) {
    if (brs.empty() || brs.front().type != cram::ReactionType::Fission)
      continue;
    const cram::Zai parent = brs.front().parent;
    if (c.chain.nearestYields(parent, 0.0253) == nullptr)
      out.push_back(parent);
  }
  return out;
}

// Nuclides whose entire decay is spontaneous fission (every mode isFission).
inline std::vector<cram::Zai> decaysOnlyBySpontaneousFission(const VeraCase& c) {
  std::vector<cram::Zai> out;
  for (const cram::Zai& z : c.chain.nuclides()) {
    const cram::DecayData* d = c.chain.decay(z);
    if (d == nullptr || d->decayConstant == 0.0 || d->modes.empty())
      continue;
    bool all = true;
    for (const cram::DecayMode& m : d->modes)
      all = all && m.isFission;
    if (all)
      out.push_back(z);
  }
  return out;
}

}  // namespace cram_validation
