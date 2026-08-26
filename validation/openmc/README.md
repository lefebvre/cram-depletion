# OpenMC reference-data generator (offline)

`generate_vera_pin.py` produces the VERA pin reference data that
[tests/validation/vera_pin_test.cpp](../../tests/validation/vera_pin_test.cpp)
replays. It is an **offline, one-time** tool: it depends on OpenMC and a
multi-GB nuclear-data library and is deliberately **not** part of the C++ build
or CI. The C++ validation test SKIPs when the generated data is absent, so a
clean checkout always builds and passes without any of this.

## Why a transport code is needed here

This repository is a transport-free depletion *engine*. The full VERA depletion
benchmark (Yu & Forget, *Ann. Nucl. Energy* 170 (2022) 108973) couples neutron
transport to the depletion solver: transport supplies the one-group cross
sections the burnup matrix is built from. We cannot generate those without a
transport code, so we let OpenMC do one beginning-of-life transport solve to
produce a **fixed** one-group micro-cross-section set, then deplete at a fixed
flux. With the burnup matrix thus held constant, the C++ engine, given the same
chain + micro xs + flux, must reproduce OpenMC's number densities to CRAM
precision — verifying our matrix assembly + CRAM against OpenMC's, isolated from
transport/cross-section error.

## Setup

OpenMC is **not** pip-installable; use conda-forge (the committed `vera_pin1a`
data was generated with OpenMC 0.15.3 installed this way):

```bash
# OpenMC + Python deps (micromamba/conda/mamba all work):
micromamba create -p ./omc -c conda-forge "openmc>=0.15" python=3.12

# Nuclear data: ENDF/B-VIII.0 HDF5 library (3.0 GB, temps 20/294/600 K):
curl -L -o endfb80.tar.xz \
  "https://zenodo.org/records/8410375/files/endfb80-lowtemp.tar.xz?download=1"
tar xf endfb80.tar.xz
export OPENMC_CROSS_SECTIONS="$PWD/endfb80-lowtemp/cross_sections.xml"

# Depletion chain: the simplified CASL PWR chain (~1 MB, 228 nuclides), as used
# by the paper. (The full 3820-nuclide chain also works but is 27 MB.) Grab it
# from the "Simplified Chain" links at https://openmc.org/data/, e.g.
curl -L -o chain_simple.xml \
  "https://anl.box.com/shared/static/3nvnasacm2b56716oh5hyndxdyauh5gs.xml"
```

`get_microxs_and_flux` spawns the `openmc` binary via PATH, so make sure the env
is activated (or `PATH`/`LD_LIBRARY_PATH` point at it) before running.

## Generate

```bash
# from the repository root (VERA_PARTICLES etc. trade statistics for speed)
VERA_PARTICLES=1500 VERA_BATCHES=25 \
  python validation/openmc/generate_vera_pin.py \
    --case 1a --chain chain_simple.xml --out tests/validation/data/vera_pin1a
```

This writes `tests/validation/data/vera_pin1a/` with `chain.xml`, `schedule.csv`,
`micro_xs.csv`, and `density.csv` (~1 MB total). The regression test then
reproduces OpenMC's number densities:

```bash
cmake --build build
ctest --test-dir build -R VeraPin --output-on-failure
```

## Notes

* The OpenMC depletion Python API changes between releases. This script targets
  the 0.15.x API (`get_microxs_and_flux`, `IndependentOperator`,
  `PredictorIntegrator`); adjust call signatures to your installed version if
  needed.
* The model is a reduced-fidelity reflective pin cell (low particle/batch counts)
  — enough for a self-consistent depletion trajectory, not for publication-grade
  statistics. Raise `VERA_PARTICLES`/`VERA_BATCHES` (or the defaults in
  `build_model()`) for tighter cross sections.
* All materials run at **294 K**: the free low-temperature library only carries
  20/294/600 K, and this is an engine-reproduction test (densities, not absolute
  VERA k-eff), so temperature fidelity is irrelevant. Consequently cases 1A and
  1C (which differ only in fuel temperature) produce identical data here — use a
  full-temperature library if you need the real per-case distinction.
* Geometry constants follow the VERA spec (pellet r=0.4096 cm, clad
  0.418/0.475 cm, pitch 1.26 cm). Enrichment/temperature per case are in `CASES`.
* Reproduction is exact by construction: in `source-rate` mode OpenMC's per-atom
  rate is `xs[b] · 1e-24 · flux` (the material volume cancels), which is exactly
  what the C++ engine assembles. Only reactions present in the chain are applied,
  and `(n,gamma)` etc. split to ground/metastable by the chain's branching ratios.
