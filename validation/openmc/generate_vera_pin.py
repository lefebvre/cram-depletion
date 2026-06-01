#!/usr/bin/env python3
"""Generate fixed-cross-section VERA pin reference data for the cram-depletion
validation suite.

This is an OFFLINE, one-time data generator -- it is NOT part of the C++ build
or CI. It needs OpenMC and a multi-GB nuclear-data library installed (see
README.md in this directory). It produces, for each VERA pin case, the reference
files consumed by tests/validation/vera_pin_test.cpp:

    tests/validation/data/<case>/chain.xml      OpenMC depletion chain
    tests/validation/data/<case>/schedule.csv   step,dt_seconds,flux
    tests/validation/data/<case>/micro_xs.csv   step,nuclide,reaction,xs_barn
    tests/validation/data/<case>/density.csv    step,nuclide,atoms

Approach (why this reproduces OpenMC in the C++ engine):
  * One transport solve at beginning-of-life produces a FIXED one-group micro
    cross-section set (openmc.deplete.MicroXS.from_model).
  * Depletion runs with IndependentOperator in 'source-rate' mode at a FIXED
    scalar flux and the first-order PredictorIntegrator. The burnup matrix is
    therefore constant and identical to what the C++ engine assembles from the
    same chain + micro xs + flux, so the C++ predictor march reproduces these
    number densities to CRAM precision. This isolates the depletion engine
    (matrix assembly + CRAM) from transport / cross-section error.

Run from the repository root, e.g.:
    python validation/openmc/generate_vera_pin.py --case 1a
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path

import numpy as np
import openmc
import openmc.deplete


# VERA depletion benchmark pin geometry (CASL-U-2015-1014-000). Common to all
# pin cases; only the fuel enrichment / temperature change between cases.
PELLET_RADIUS = 0.4096  # cm
CLAD_INNER = 0.418      # cm
CLAD_OUTER = 0.475      # cm
PIN_PITCH = 1.26        # cm
FUEL_DENSITY = 10.257   # g/cm^3 (UO2)

# Reactions we export (matching the C++ ReactionType set). OpenMC reaction
# names map directly to the strings the C++ reader expects.
REACTIONS = ["fission", "(n,gamma)", "(n,2n)", "(n,3n)", "(n,4n)", "(n,a)", "(n,p)"]

CASES = {
    "1a": dict(enrichment=3.1, fuel_temp=565.0),
    "1c": dict(enrichment=3.1, fuel_temp=900.0),
}


def build_model(enrichment: float, fuel_temp: float) -> openmc.Model:
    """A single reflective PWR fuel pin cell."""
    # All materials at 294 K. The benchmark's per-case fuel temperature is not
    # used: the freely-available low-temperature HDF5 library only carries
    # 20/294/600 K, and this is an engine-reproduction test (the depletion
    # densities, not absolute VERA k-eff, are what we compare), so temperature
    # fidelity is irrelevant here.
    temp = 294.0
    fuel = openmc.Material(name="fuel")
    fuel.add_element("U", 1.0, enrichment=enrichment)
    fuel.add_element("O", 2.0)
    fuel.set_density("g/cm3", FUEL_DENSITY)
    fuel.temperature = temp
    fuel.volume = np.pi * PELLET_RADIUS**2  # per unit height; needed for depletion
    fuel.depletable = True

    clad = openmc.Material(name="clad")
    clad.add_element("Zr", 1.0)
    clad.set_density("g/cm3", 6.55)
    clad.temperature = temp

    water = openmc.Material(name="water")
    water.add_element("H", 2.0)
    water.add_element("O", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")
    water.set_density("g/cm3", 0.743)
    water.temperature = temp

    materials = openmc.Materials([fuel, clad, water])

    r_fuel = openmc.ZCylinder(r=PELLET_RADIUS)
    r_clad_in = openmc.ZCylinder(r=CLAD_INNER)
    r_clad_out = openmc.ZCylinder(r=CLAD_OUTER)
    half = PIN_PITCH / 2.0
    box = openmc.model.RectangularPrism(PIN_PITCH, PIN_PITCH, boundary_type="reflective")

    fuel_cell = openmc.Cell(fill=fuel, region=-r_fuel)
    gap_cell = openmc.Cell(region=+r_fuel & -r_clad_in)  # void gap
    clad_cell = openmc.Cell(fill=clad, region=+r_clad_in & -r_clad_out)
    water_cell = openmc.Cell(fill=water, region=+r_clad_out & -box)
    geometry = openmc.Geometry([fuel_cell, gap_cell, clad_cell, water_cell])

    settings = openmc.Settings()
    settings.particles = int(os.environ.get("VERA_PARTICLES", 5000))
    settings.batches = int(os.environ.get("VERA_BATCHES", 40))
    settings.inactive = int(os.environ.get("VERA_INACTIVE", 10))
    settings.run_mode = "eigenvalue"

    return openmc.Model(geometry, materials, settings)


def export_case(case: str, out_dir: Path, flux: float, days, chain_file: str) -> None:
    params = CASES[case]
    model = build_model(params["enrichment"], params["fuel_temp"])
    fuel = next(m for m in model.materials if m.name == "fuel")

    # One BOL transport solve -> a FIXED one-group micro cross-section set.
    # A single energy group (0..20 MeV) collapses the spectrum so the exported
    # xs is a plain barn value: reaction rate = xs[b] * 1e-24 * flux[n/cm^2/s],
    # the same relation the C++ engine uses.
    fluxes, micros = openmc.deplete.get_microxs_and_flux(
        model, [fuel], reactions=REACTIONS, energies=[0.0, 20.0e6], chain_file=chain_file
    )
    micro_xs = micros[0]  # MicroXS: data[nuclide, reaction, group], one group

    # Fixed-flux ('source-rate') depletion with the first-order predictor so the
    # burnup matrix is constant and reproducible by the C++ engine. The flux we
    # hand the operator is volume-integrated [n-cm/s]; source_rate carries the
    # magnitude so the per-atom rate is xs * flux.
    operator = openmc.deplete.IndependentOperator(
        openmc.Materials([fuel]),
        fluxes=[np.array([flux * fuel.volume])],
        micros=[micro_xs],
        chain_file=chain_file,
        normalization_mode="source-rate",
    )
    integrator = openmc.deplete.PredictorIntegrator(
        operator, days, source_rates=1.0, timestep_units="d"
    )
    integrator.integrate()
    results = openmc.deplete.Results("depletion_results.h5")

    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(chain_file, out_dir / "chain.xml")

    dt_seconds = [d * 86400.0 for d in days]
    with open(out_dir / "schedule.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "dt_seconds", "flux"])
        for k, dt in enumerate(dt_seconds, start=1):
            w.writerow([k, dt, flux])

    # Fixed micro xs (data[i, j, 0]); constant over the trajectory, so written
    # once. The C++ reader applies it to every step.
    nuclides = list(micro_xs.nuclides)
    reactions = list(micro_xs.reactions)
    with open(out_dir / "micro_xs.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nuclide", "reaction", "xs_barn"])
        for i, nuc in enumerate(nuclides):
            for j, rxn in enumerate(reactions):
                xs = float(micro_xs.data[i, j, 0])
                if xs > 0.0:
                    w.writerow([nuc, rxn, repr(xs)])

    # Reference number densities (atoms) per step, every depleted nuclide.
    with open(out_dir / "density.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "nuclide", "atoms"])
        for nuc in nuclides:
            try:
                _, atoms = results.get_atoms(fuel, nuc)
            except KeyError:
                continue
            for step, n in enumerate(atoms):
                w.writerow([step, nuc, repr(float(n))])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", default="1a", choices=sorted(CASES), help="VERA pin case")
    ap.add_argument("--flux", type=float, default=3.0e14, help="fixed scalar flux [n/cm^2/s]")
    ap.add_argument(
        "--chain", default=os.environ.get("OPENMC_CHAIN_FILE", "chain.xml"),
        help="OpenMC depletion chain XML (or set OPENMC_CHAIN_FILE / openmc.config)",
    )
    ap.add_argument(
        "--out", default=None,
        help="output directory (default tests/validation/data/vera_pin<case>)",
    )
    args = ap.parse_args()

    # VERA-style coarse '20-step' schedule (EFPD): 0.25, 6.0, 6.25, 12.5,
    # 2x25, 14x100. Truncated here -- extend to the full 1500 EFPD if desired.
    days = [0.25, 6.0, 6.25, 12.5, 25.0, 25.0] + [100.0] * 8

    out_dir = Path(args.out) if args.out else Path(
        f"tests/validation/data/vera_pin{args.case}"
    )
    export_case(args.case, out_dir, args.flux, days, args.chain)
    print(f"wrote reference data to {out_dir}")


if __name__ == "__main__":
    main()
