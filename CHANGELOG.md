# Changelog

Releases are tagged `vMAJOR.MINOR.PATCH`; the version lives in `project()` in
the top-level `CMakeLists.txt` and CI refuses a tag that disagrees with it.
Past 1.0 the major is the breaking-change axis, so a consumer that asked for
`find_package(cram 2.0)` is satisfied by any 2.x.

## 2.1.0

Burnup: everything upstream of the matrix that a fixed-cross-section depletion
needs, and the adjoint of the linear problem.

- `cram/integrator.hpp`: Predictor, CE/CM, CE/LI, LE/QI and CF4 integrators
  marching `dn/dt = A(n) n` over a step schedule, with coefficients and
  application order matching `openmc.deplete`.
- `cram/deplete.hpp`: `DepletionSystem`, one material with fixed one-group
  cross sections under constant flux or constant power. Refuses inputs that
  would fail silently: a fission channel with no energy release under a power
  target, a negative cross section, a composition with no fissile material
  asked to hold a non-zero power.
- `cram/reaction.hpp`: reaction channels, their ground-state products, the
  OpenMC reaction names, and `ChainReaction` for a chain's reaction topology.
- `cram/chain_xml.hpp` (optional, `CRAM_WITH_CHAIN_XML`, pugixml): reader for
  OpenMC `depletion_chain` XML. Unusable entries are dropped and counted, never
  stored as phantom nuclides.
- `cram/adjoint.hpp`: adjoint solves and marches for linear systems,
  first-order sensitivities `dR/dλ` and `dR/dσ` by quadrature over the
  trajectories, time-integrated responses, and the generators for
  `dn/dt = A n + s` and for exact inventory integrals.
- `DecayMode::daughter`: a decay mode may name its product explicitly; matrix
  assembly and `close()` honor it. `DecayMode` initializations should be
  designated so the new trailing member is omitted deliberately.
- `tests/validation/`: replay of an OpenMC VERA pin depletion; the engine
  reproduces OpenMC's number densities to below 1e-3 for the benchmark
  nuclides given identical chain, cross sections and flux.

Additive: 2.0 consumers build unchanged.

## 2.0.0

Source-compatibility break: the data structs no longer default-initialize
members whose zero value would be a plausible wrong answer. `DecayMode`
requires all of `rtyp`, `branching`, `finalState` and `isFission`;
`DecayData::halfLife`, `FissionYields::energy` and `Zai::z`/`Zai::a` must be
supplied at every braced initialization. Solver results are unchanged.

- Chain ordering made independent of hash iteration order, so the matrix
  layout and the last bits of every result are a function of the data alone.
- `CramSolver` reuses the symbolic analysis across `prepare()` calls with the
  same sparsity pattern, and exposes `beginPrepare()/preparePole()/endPrepare()`
  so a caller can drive the per-pole factorizations in parallel.
- Per-region matrix assembly cost cut; benchmarks cover the decay-only and
  burnup regimes.
- Install and export configuration; the version centralized on
  `project(VERSION)` with `cram/version.hpp` and `deplete --version`.
- Every target compiled with warnings on and CI gated on them.
- `CRAM_WITH_ENDFTK=ON` fixed against consteval format-string checks.

## 1.0.1

- MSVC build fixed and a Windows CI job added.
- Project CMake options prefixed `CRAM_`; Eigen's own test suite kept out of
  CTest.
- GoogleTest discovery deferred to test time for multi-config reliability.

## 1.0.0

Initial release: IPF-CRAM16/48 solver, `CramSolver` with cached per-pole
factorizations, `DepletionChain` with decay data, fission yields and matrix
assembly, and the optional ENDFtk reader.
