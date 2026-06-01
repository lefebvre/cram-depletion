# cram-depletion

[![CI](https://github.com/lefebvre/cram-depletion/actions/workflows/ci.yml/badge.svg)](https://github.com/lefebvre/cram-depletion/actions/workflows/ci.yml)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)]()
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-yellow.svg)](LICENSE.md)
[![codecov](https://codecov.io/gh/lefebvre/cram-depletion/graph/badge.svg)](https://codecov.io/gh/lefebvre/cram-depletion)

A small C++20 CRAM (Chebyshev Rational Approximation Method) solver for the
nuclear depletion / Bateman equations, with an optional ENDF/B-VIII data reader
built on [njoy/ENDFtk](https://github.com/njoy/ENDFtk).

Solves `n(t) = exp(A t) n0` where `A` is the burnup matrix
(`A(j,i)` = rate at which nuclide `i` produces nuclide `j`).

## Layout

```
cram/nuclide.hpp         ZAI identity + ENDF decay-mode -> daughter logic
cram/chain.hpp           DepletionChain: nuclides, decay data, fission yields
cram/cram.hpp            CRAM solver interface
cram/cram.cpp            IPF-CRAM-16 / CRAM-48 (verified coefficients)
cram/burnup_matrix.cpp   matrix assembly (decay, fission source, reactions)
cram/integrator.hpp      time integrators (predictor, CE/CM, CE/LI, LE/QI, CF4)
cram/deplete.hpp         constant-power depletion system (fixed one-group XS)
cram/chain_xml.cpp       OpenMC depletion_chain XML reader (via pugixml)
cram/endf_reader.cpp     ENDFtk ingestion (optional, see notes)
cram-apps/deplete.cpp    runnable demo + ENDF driver
cmake/gcov_to_lcov.py    gcov JSON -> LCOV tracefile (consumed by VS Code)
cmake/ResetCoverage.cmake  clear stale .gcda artifacts before a capture
external/                FetchContent cache (gitignored; survives clean builds)
```

## Build

Requires CMake >= 3.28 and a C++20 compiler. Eigen, googletest, and ENDFtk
are pulled via FetchContent if not found system-wide; their sources land in
`external/` so they're not re-cloned on every clean build.

```bash
# core + demo
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/deplete            # self-contained decay demo

# with the ENDFtk reader (pulls ENDFtk v1.2.0 + transitive deps)
cmake -B build -DWITH_ENDFTK=ON
cmake --build build
./build/deplete decay_sublibrary.endf 86400
ctest --test-dir build -R EndfReaderIntegration   # parses real ENDF data
```

## Tests, coverage, and quality checks

Unit tests use GoogleTest (found via `find_package`, else fetched). The CRAM
solver is checked against independent closed-form Bateman solutions, not a
second numerical method, so a solver bug cannot hide behind a matching
reference bug.

```bash
# build + run the unit + integration suite
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure

# coverage -> build/coverage/coverage.info  (LCOV tracefile)
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
cmake --build build
cmake --build build --target coverage

# Address + UB sanitizers
cmake -B build -DENABLE_SANITIZERS=ON && cmake --build build && ctest --test-dir build

# format + static analysis
clang-format --dry-run -Werror cram/*.{hpp,cpp} cram-apps/*.cpp tests/**/*.{hpp,cpp}
cmake --build build --target cppcheck      # cppcheck
clang-tidy -p build cram/*.cpp             # clang-tidy (config in .clang-tidy)
```

Coverage is produced via `cmake/gcov_to_lcov.py` (a thin wrapper around
`gcov -j -b -c`), so no `gcovr`/`lcov` install is needed locally. The output
file path matches `cmake.coverageInfoFiles` in `.vscode/settings.json`, so
the CMake Tools extension surfaces line/branch gutters automatically once
the `coverage` target has been built.

`.github/workflows/ci.yml` runs format-check (gating the rest of the
pipeline so noncompliant formatting fails fast), then build-test-coverage
(with `gcovr` enforcing a 95% line-coverage floor + HTML artifact),
sanitizers, static-analysis, and the ENDFtk integration job. Test layout:

```
tests/bateman.hpp        analytic Bateman references
tests/test_nuclide.cpp   ZAI packing, RTYP decode, decay-mode -> daughter
tests/test_chain.cpp     registration, decay/yield data, matrix assembly
tests/test_cram.cpp      CRAM16/CRAM48 vs analytic; mass, stiffness, edge cases
tests/test_integrator.cpp  integrator order-of-accuracy + constant-A exactness
tests/test_deplete.cpp     constant-power flux normalization + assembly
tests/test_chain_xml.cpp   OpenMC depletion_chain parsing
tests/validation/        guarded replay vs OpenMC-generated VERA pin data
tests/integration/       real ENDFtk reader vs real ENDF data (WITH_ENDFTK only)
```

## Validation against OpenMC (VERA depletion benchmark)

The depletion engine mirrors the validation approach of OpenMC's depletion
module (Yu & Forget, *Ann. Nucl. Energy* 170 (2022) 108973, which verifies
OpenMC against the VERA depletion benchmark). Two layers, both transport-free:

1. **Integrator order-of-accuracy** (`tests/test_integrator.cpp`, always in CI).
   The predictor / CE-CM / CE-LI / LE-QI / CF4 integrators (`cram/integrator.hpp`,
   coefficients matching OpenMC's `openmc.deplete`) are marched over a
   constant-power pin problem at refining time steps; the observed convergence
   orders (≈1 predictor, ≈2 CE/CM·CE/LI·LE/QI, ≈4 CF4) reproduce the paper's
   §4.5 / Figs 15–18. When the flux is held fixed (`A` constant) every scheme
   collapses to the exact `exp(A·dt)`, which is also checked.

2. **Engine reproduction of OpenMC** (`tests/validation/`). OpenMC runs a
   fixed-cross-section VERA pin depletion (ENDF/B-VIII.0, simplified CASL chain)
   with its predictor; given the same chain + one-group micro cross sections +
   flux, this engine's predictor march reproduces OpenMC's number densities to
   **< 1e-3** for every benchmark nuclide of interest (U/Np/Pu/Am isotopes,
   Xe-135, Cs-137, Nd-148, Sm-149, Gd-157 — the major ones to ~1e-5). This
   checks matrix assembly + CRAM, isolated from the transport / cross-section
   error a transport-free engine cannot reproduce. The committed reference data
   for case `vera_pin1a` lives under `tests/validation/data/`; regenerate or add
   cases offline with [validation/openmc/generate_vera_pin.py](validation/openmc/).
   The test SKIPs (never fails) when the data directory is absent, so a clean
   checkout without it still passes.

## How CRAM works here

`exp(At) n0 ≈ α₀ · Πₗ ( I + 2·Re( αₗ (At − θₗ I)⁻¹ ) ) n0`, applied
incrementally to the running vector (Pusa's IPF form, NSE 182:3, 2016). Each
pole is one complex sparse solve via `Eigen::SparseLU`; CRAM16 = 8 poles,
CRAM48 = 24. The θ/α/α₀ constants in `cram/cram.cpp` are transcribed from
OpenMC (MIT license) and were checked against analytic Bateman solutions to
~1e-15 relative error.

## Things that will bite you (physics / data)

- **Independent vs cumulative yields.** Use MF8/MT454 (independent) yields when
  you model the decay chain explicitly. MT459 (cumulative) yields already fold
  in decay feeding — combining them with an explicit chain double-counts. The
  reader defaults to MT454.
- **Fission only contributes if there is a fission rate.** Yields are a
  branching distribution, not a source by themselves. Spontaneous fission is
  driven by the SF decay branch (RTYP 6) + spontaneous-fission yields (energy 0);
  neutron-induced fission needs `fissionRate = σ_f · φ` fed via
  `addFissionSource(...)`.
- **Register daughters.** A production term is dropped if the daughter isn't in
  the chain. `loadDecayData` adds reachable daughters; if you build chains by
  hand, `chain.add(...)` every product first.
- **CRAM accuracy is absolute, not relative, for trace species.** Nuclides many
  orders of magnitude below the dominant one can come out slightly negative
  (~1e-17). That is expected; clamp at zero if you need non-negativity.
- **ENDFtk API.** The reader targets ENDFtk's current tree + parsed-section
  interface. The accessor names are verified against the ENDFtk source and a
  real build: `tests/integration` parses a real MT457 decay section (Am-242m)
  and a real MT454 fission-yield file (Pu-239) and checks the results. Notable
  API facts pinned down there: `halfLife()` returns `{value, uncertainty}`;
  `branchingRatio()` returns a `{value, uncertainty}` range; the parent
  isomeric state is `LISO()`; fission products come from
  `fissionProductIdentifiers()` / `isomericStates()` / `fissionYieldValues()`.
  MF8 supports MT457 (decay) and MT454/459 (fission yields); older docs that
  say "MT457 only" are out of date.

## Possible next steps

- The stochastic-implicit (SI-CE/LI, SI-LE/QI) and EPC-RK4 integrators from the
  OpenMC set, for Xe-stability on very large systems.
- One-group collapse of MF3 cross sections (or ACE data) to get reaction rates
  for `(n,γ)`, `(n,2n)`, fission, etc.
- A sparsity-preserving ordering / reuse of the symbolic factorization across
  poles for speed on the full ~3800-nuclide ENDF/B-VIII chain.
