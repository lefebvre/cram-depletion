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
cram/endf_reader.cpp     ENDFtk ingestion (optional, see notes)
cram-apps/deplete.cpp    runnable demo + ENDF driver
cmake/gcov_to_lcov.py    gcov JSON -> LCOV tracefile (consumed by VS Code)
cmake/ResetCoverage.cmake  clear stale .gcda artifacts before a capture
external/                FetchContent cache (gitignored; survives clean builds)
```

## Build

Requires CMake >= 3.28 and a C++20 compiler. Eigen, googletest, and ENDFtk
are pulled via FetchContent if not found system-wide; their sources land in
`external/` so they're not re-cloned on every clean build. The Eigen fallback
defaults to v5.0.1 and is configurable (see `CRAM_EIGEN_GIT_TAG` below).

```bash
# core + demo
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/deplete            # self-contained decay demo

# with the ENDFtk reader (pulls ENDFtk v1.2.0 + transitive deps)
cmake -B build -DCRAM_WITH_ENDFTK=ON
cmake --build build
./build/deplete decay_sublibrary.endf 86400
ctest --test-dir build -R EndfReaderIntegration   # parses real ENDF data
```

All project build options are `CRAM_`-prefixed (`CRAM_WITH_ENDFTK`,
`CRAM_ENABLE_TESTS`, `CRAM_ENABLE_COVERAGE`, `CRAM_ENABLE_SANITIZERS`,
`CRAM_ENABLE_TSAN`, `CRAM_ENABLE_BENCHMARKS`, `CRAM_ENABLE_CLANG_TIDY`,
`CRAM_WERROR`, `CRAM_ENABLE_INSTALL`) so they don't collide with the options of
dependencies fetched via FetchContent.

Which Eigen gets used is controlled by three cache variables:

| Variable | Default | Meaning |
| --- | --- | --- |
| `CRAM_EIGEN_MIN_VERSION` | `3.4` | A system Eigen older than this is ignored and the fetch fallback is used instead |
| `CRAM_EIGEN_GIT_REPOSITORY` | `https://gitlab.com/libeigen/eigen.git` | Repository cloned when Eigen must be fetched |
| `CRAM_EIGEN_GIT_TAG` | `5.0.1` | Tag, branch, or commit fetched. A bare commit SHA disables the shallow clone automatically |

`find_package(Eigen3)` is deliberately called *without* a version argument:
Eigen's package version file uses `SameMajorVersion` compatibility, so asking
for 3.4 would reject a system Eigen 5.x. The floor is applied afterwards from
`CRAM_EIGEN_MIN_VERSION`; the configure log always reports which Eigen won.

```bash
# build against a different Eigen (ignores any system install)
cmake -B build -DCRAM_EIGEN_GIT_TAG=3.4.0 -DCMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON
```

To exercise the fetch path when a system Eigen is present, use CMake's stock
escapes — `-DCMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON` or
`-DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER`. No project-specific option is
needed for that.

## Installing and consuming

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/cram
cmake --build build
cmake --install build
```

That installs the `cram` static library, the public headers under
`include/cram/`, the `deplete` driver, and a package config, after which a
downstream project needs only:

```cmake
find_package(cram 1.0 REQUIRED)
target_link_libraries(my_target PRIVATE cram::cram)
```

`#include <cram/cram.hpp>` then resolves from the install tree.
`cram_poles_internal.hpp` is not installed: it is an implementation detail
included only by `cram/*.cpp`.

### Versioning

The version lives in exactly one place — the `project()` call in the top-level
`CMakeLists.txt`. Everything else derives from it:

| Consumer | How it gets the version |
| --- | --- |
| `cram/version.hpp` | Generated into the build tree from `cram/version.hpp.in`; provides `CRAM_VERSION_{MAJOR,MINOR,PATCH,STRING,HEX}` and `cram::kVersion` |
| `find_package(cram 1.0)` | `cramConfigVersion.cmake`, generated with `SameMajorVersion` compatibility |
| `deplete --version` | Prints `cram::kVersion` |
| release tag | The `version-tag-check` CI job fails a `v*` tag that disagrees with `project()` |

So cutting a release is: bump `project(cram_depletion VERSION ...)`, commit,
then tag to match. Tagging without the bump fails CI rather than shipping a
binary whose `--version` contradicts its tag.

```cpp
#include <cram/version.hpp>
#if CRAM_VERSION_HEX < 0x010000
#error "needs cram >= 1.0.0"
#endif
```

**Installing requires a system Eigen.** `CRAM_ENABLE_INSTALL` defaults to `ON`
only for a top-level build that found Eigen via `find_package`. A *fetched*
Eigen is a target in the build tree that is never itself installed, so an
exported `cram` naming it cannot produce a working package — CMake rejects it
outright with `install(EXPORT "cramTargets" ...) includes target "cram" which
requires target "eigen" that is not in any export set`. Rather than fail the
configure for the many developers who build without a system Eigen and never
install, the option simply defaults off there; setting it `ON` anyway reports
the reason and the remedy. The same applies to `CRAM_WITH_ENDFTK`, since ENDFtk
is always fetched, so the two options are mutually exclusive for now.

`CRAM_ENABLE_INSTALL` also defaults off whenever cram is consumed via
`add_subdirectory` or `FetchContent`, so cram's headers and package config
don't leak into the parent project's install tree.

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
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCRAM_ENABLE_COVERAGE=ON
cmake --build build
cmake --build build --target coverage

# Address + UB sanitizers
cmake -B build -DCRAM_ENABLE_SANITIZERS=ON && cmake --build build && ctest --test-dir build

# format + static analysis
clang-format --dry-run -Werror cram/*.{hpp,cpp} cram-apps/*.cpp tests/**/*.{hpp,cpp}
cmake --build build --target cppcheck      # cppcheck
clang-tidy -p build cram/*.cpp             # clang-tidy (config in .clang-tidy)
```

CTest test names are prefixed with `cram.` (e.g. `cram.Cram.ScalarExponential`)
to namespace this project's tests apart from any dependency's. Eigen 3.4, when
pulled via FetchContent rather than found system-wide, would otherwise register
its own ~900-test suite into the same CTest project; `BUILD_TESTING` is forced
off around Eigen's configuration in `CMakeLists.txt` (and the user's own value
restored afterwards) so `ctest` only sees the CRAM suite. Eigen >= 5 gates its
suite on `PROJECT_IS_TOP_LEVEL` and needs no such help, but 3.4-era tags remain
selectable via `CRAM_EIGEN_GIT_TAG`, so the guard stays.

Coverage is produced via `cmake/gcov_to_lcov.py` (a thin wrapper around
`gcov -j -b -c`), so no `gcovr`/`lcov` install is needed locally. The output
file path matches `cmake.coverageInfoFiles` in `.vscode/settings.json`, so
the CMake Tools extension surfaces line/branch gutters automatically once
the `coverage` target has been built.

`.github/workflows/ci.yml` runs format-check (gating the rest of the
pipeline so noncompliant formatting fails fast), then build-test-coverage
(with `gcovr` enforcing a 95% line-coverage floor + HTML artifact),
a Windows MSVC build-and-test job, sanitizers, static-analysis, and the
ENDFtk integration job. Test layout:

```
tests/bateman.hpp        analytic Bateman references
tests/test_nuclide.cpp   ZAI packing, RTYP decode, decay-mode -> daughter
tests/test_chain.cpp     registration, decay/yield data, matrix assembly
tests/test_cram.cpp      CRAM16/CRAM48 vs analytic; mass, stiffness, edge cases
tests/integration/       real ENDFtk reader vs real ENDF data (WITH_ENDFTK only)
```

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
- **Build the data structs with braces.** `Zai`, `DecayMode`, `DecayData` and
  `FissionYields` leave the members whose zero value would be a plausible wrong
  answer — `Zai::z`/`::a`, every `DecayMode` field, `DecayData::halfLife`,
  `FissionYields::energy` — without a default initializer, so omitting one is a
  `-Wmissing-field-initializers` diagnostic rather than a nuclide that silently
  never decays or a yield table that silently becomes a spontaneous-fission one.
  `DecayData d; d.halfLife = ...;` leaves the rest indeterminate: write
  `DecayData d{.halfLife = ...}` instead. Members that do keep an initializer —
  `Zai::i` (ground state), `DecayData::decayConstant` (derived by `setDecay`),
  `gammaEnergyPerDecay`, and the two vectors — are optional by design.
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

- Predictor-corrector time integration (CE/CM, LE/QI) for coupling to a flux
  solver, so the matrix is rebuilt as the spectrum/number densities change.
- One-group collapse of MF3 cross sections (or ACE data) to get reaction rates
  for `(n,γ)`, `(n,2n)`, fission, etc.
- A sparsity-preserving ordering / reuse of the symbolic factorization across
  poles for speed on the full ~3800-nuclide ENDF/B-VIII chain.
