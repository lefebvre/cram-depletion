# Changelog

Releases are tagged `vMAJOR.MINOR.PATCH`; the version lives in `project()` in
the top-level `CMakeLists.txt` and CI refuses a tag that disagrees with it.
Past 1.0 the major is the breaking-change axis, so a consumer that asked for
`find_package(cram 2.0)` is satisfied by any 2.x. "Breaking" is judged against
cram's actual consumers, which today means nufall: a change none of them can
observe ships in a minor release even when it tightens behaviour, and the entry
marks it **stricter** so that anyone else relying on the old behaviour can see
it.

## 2.3.0

Review fixes: two data bugs, stricter input handling, block-triangular
factorization, and a validation suite held to round-off.

Every function keeps its signature. Code that copies or moves a
`DepletionSystem`, or tests the old `WITH_ENDFTK` macro, no longer compiles, and
input that 2.2 accepted and quietly mishandled now throws; those entries are
marked **stricter**. nufall uses none of the affected paths, so this is a minor
release.

The `cramPoles()` coefficients published in 2.2.0 are unchanged. 2.2.0's
"results are bit-identical" no longer holds for `cramSolve()` and `CramSolver`:
the new factorization order moves their results at round-off level.

### Fixes

- ENDF fission yields of an isomeric parent (e.g. Am-242m1) were filed under
  its ground state, because MT454/459 carry no LISO. The parent's state now
  comes from MF1/MT451, and a material without that header is refused.
- Spontaneous and neutron-induced yields no longer stand in for each other.
  Spontaneous fission uses only the table at energy 0. `nearestYields()` at a
  non-zero energy considers only neutron-induced tables, and at energy 0
  returns the SF table or null, never the nearest thermal table. An SF mode
  whose parent has no SF table removes the parent as before and produces
  nothing; that loss is now reported, as an unregistered daughter is.
- `decayConstantSensitivities()` read S with no checks, so an S built for
  another chain was an out-of-bounds read, and an entry S never integrated
  read as zero. It now rejects an S_k of the wrong size, or one missing an
  entry the decay needs, checking each interval separately so that one cannot
  hide another's gap. `reactionSensitivities()` already made both checks.
- The FetchContent base directory was shared by every build tree, so trees
  overwrote each other's dependency builds. Dependencies now fetch into each
  build tree. Paths are project-scoped, so `add_subdirectory` consumers work.

### Stricter

- Built without ENDFtk, `loadDecayData` and `loadFissionYields` throw
  `std::runtime_error` instead of returning 0, as the chain XML loaders already
  did. The compile definition is renamed from `WITH_ENDFTK` to
  `CRAM_WITH_ENDFTK`.
- The chain XML reader throws on a malformed number (including a doubled sign
  such as `+-0.5`), a `<fission_yields>` table without an energy, or one whose
  products and data differ in length, instead of guessing. Numbers are parsed
  independently of the C locale. `parseNuclideName` accepts only the whole
  `<symbol><mass>[_m<level>]` grammar.
- `DepletionChain::add` throws for a `Zai` its packed key cannot hold, and
  `Zai::operator==` compares fields rather than keys, so `{1, 2, 10}` and
  `{1, 3, 0}` no longer compare equal.
- `cramSolve`, `CramSolver::prepare` and the integrators throw
  `std::invalid_argument` for a negative or non-finite time step.
- `DepletionSystem` is neither copyable nor movable, since `matrixBuilder()`
  binds its address. `intervalMatrices()` now restores the system's flux on
  return.

### Added

- `requireTimeStep(dt)` in `cram/cram.hpp`: the check every solver entry point
  applies.
- `Zai::representable()`: whether a `Zai` fits its packed key
  (`z >= 0`, `1 <= a <= 999`, `0 <= i <= 9`).
- `DepletionChain::spontaneousFissionYields()`, and a second, defaulted
  out-parameter of `decayMatrix()` counting SF modes with no SF table; left
  null, the count is a warning on stderr.
- `DepletionChain::decayStencil()` and `reactionStencil()`: the matrix entries
  of each process, stated once. Matrix assembly, the interval-matrix structure
  and both sensitivity contractions are built on them.
- `CramSolver::applyTranspose()`: exp(A^T dt) from the factorizations a
  `prepare(A, dt)` already holds.

### Performance

- The pole matrices are permuted to block-triangular form, with COLAMD inside
  each block, in both `cramSolve()` and `CramSolver`. At N = 1675 with CRAM48:
  `CramSolver::prepare` 31 -> 16 ms on the burnup shape and 19.5 -> 9.3 ms on
  the decay-only shape; L+U fill 2.24x -> 1.88x and 1.70x -> 1.40x (no fill-in
  at all for decay-only); a 100-step burnup march 1870 -> 1386 ms. `apply()`
  is ~10% slower on the burnup shape and ~10% faster on decay-only.
- `rateSensitivities()` serves a forward quadrature node and its mirror
  adjoint node from one factorization, through `applyTranspose()`, halving its
  factorizations: 825 -> ~350 ms on the burnup shape at N = 1024.

### Validation, tests and CI

- The OpenMC VERA replay compares every reported nuclide at every step against
  a relative 1e-10 plus an absolute 1e-12 of the largest inventory (it was
  1e-3, at end of life only), and holds nuclides the reference reports as zero
  to that floor.
- `validation/openmc/generate_constant_power_refs.py` and the
  `VeraPinConstantPower` tests compare every integrator against OpenMC at
  constant power. No reference data is committed yet, so those tests skip.
- Compile-fail tests check that omitting a required struct field is diagnosed
  (GCC and Clang; MSVC has no such warning).
- clang-tidy's bugprone and performance findings fail the static-analysis job.
  GitHub Actions are pinned to commits.

## 2.2.0

Public CRAM coefficients, for callers that bring their own linear algebra.

- `cram/cram_poles.hpp`: the CRAM16 and CRAM48 coefficients are public, through
  `cramPoles(CramOrder)`, for callers that evaluate the same IPF rational
  approximation with their own linear algebra, for example a substitution on a
  decay-only matrix, which is triangular in decay order. The header includes no
  Eigen header. It replaces the uninstalled `cram_poles_internal.hpp`, and
  `cramSolve()` and `CramSolver` now read their tables through `cramPoles()`, so
  a caller cannot drift from the library's own coefficients. The tables are
  unchanged byte for byte, and results are bit-identical.
- `CramOrder` moves to `cram/cram_poles.hpp`. `cram/cram.hpp` includes it, so
  existing code that names `cram::CramOrder` compiles unchanged.

## 2.1.0

Burnup: everything upstream of the matrix that a fixed-cross-section depletion
needs, and the adjoint of the linear problem.

- `cram/integrator.hpp`: Predictor, CE/CM, CE/LI, LE/QI and CF4 integrators
  marching `dn/dt = A(n) n` over a step schedule, with coefficients and
  application order matching `openmc.deplete`. Every exponential of a march
  goes through one `CramSolver` held by the integrator, so the symbolic
  analysis of the pole factorizations is paid once per march rather than once
  per exponential (27% off each one at N=1675); an integrator therefore holds K
  factorizations while it lives.
- `cram/deplete.hpp`: `DepletionSystem`, one material with fixed one-group
  cross sections under constant flux or constant power. Refuses inputs that
  would fail silently: a fission channel with no energy release under a power
  target, a negative cross section, a negative or non-finite flux or power, a
  composition with no fissile material asked to hold a non-zero power, and a
  chain modified after the system cached its decay matrix (`refreshChain()`
  re-caches it). A fissionable nuclide with no yield table is consumed by its
  fission channel like any other reaction with an untracked product, so the
  power normalization cannot pay for fissions the matrix never performs.
- `cram/reaction.hpp`: reaction channels, their ground-state products, the
  OpenMC reaction names, and `ChainReaction` for a chain's reaction topology.
- `cram/chain_xml.hpp` (optional, `CRAM_WITH_CHAIN_XML`, pugixml): reader for
  OpenMC `depletion_chain` XML, including the `<neutron_fission_yields
  parent="..."/>` delegation form and the target OpenMC names on a `<decay
  type="sf">` mode. A reaction target the file names but never declares is
  registered as a bare terminator, as `close()` does for a decay daughter, so
  no production is silently dropped. Unusable entries are dropped and counted,
  never stored as phantom nuclides; an absent or `"Nothing"` target is the
  file's own statement that the product is untracked, and is not counted as a
  failure.
- `cram/adjoint.hpp`: adjoint solves and marches for linear systems,
  first-order sensitivities `dR/dλ` and `dR/dσ` by quadrature over the
  trajectories, time-integrated responses, and the generators for
  `dn/dt = A n + s` and for exact inventory integrals. `intervalMatrices()`
  carries a structural entry for every channel of the system, including one
  with a zero cross section, so its derivative is integrated rather than
  half-formed; `reactionSensitivities()` refuses an `S` that lacks one.
- `DecayMode::daughter`: a decay mode may name its product explicitly; matrix
  assembly and `close()` honor it. `DecayMode` initializations should be
  designated so the new trailing member is omitted deliberately.
- `tests/validation/`: replay of an OpenMC VERA pin depletion. Given identical
  chain, cross sections and flux, the engine reproduces OpenMC's number
  densities to below 1e-3 for every nuclide the reference reports as present --
  all 181 of case `vera_pin1a`, asserted one by one, worst 4e-13 — not a
  curated subset.
- `validation/report/`: a `validation-report` build target that replays the case
  through the engine and writes the whole comparison — every nuclide, every step
  — as JSON, and a standard-library-only Python renderer that turns that into a
  standalone HTML report with inline SVG figures. Measurement and rendering are
  separate steps, so the build's only hard dependency stays on the C++ side and
  the JSON remains the machine-readable record. CI regenerates both on every
  push and publishes them as an artifact; `vera_report --fail-over-tolerance N`
  turns the comparison into a gate.

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
