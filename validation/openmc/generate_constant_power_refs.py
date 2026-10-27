#!/usr/bin/env python3
"""Generate constant-power reference trajectories for every depletion integrator.

This is an OFFLINE data generator -- it is NOT part of the C++ build or CI. It
needs OpenMC's Python API (see README.md in this directory) but, unlike
generate_vera_pin.py, no transport run and no nuclear-data library: it rebuilds
the depletion problem from a case that generator already wrote, so the inputs
are exactly the ones the C++ test reads.

    tests/validation/data/<case>/constant_power/power.csv     power_w
    tests/validation/data/<case>/constant_power/<scheme>.csv  step,nuclide,atoms

for <scheme> in predictor, cecm, celi, leqi, cf4. They are consumed by the
VeraPinConstantPower test in tests/validation/vera_pin_test.cpp, which SKIPs
while they are absent.

Why this case and not the existing one: the committed vera_pin1a trajectory is
a fixed flux with fixed cross sections, so its burnup matrix never changes and
every integrator collapses to the same exact exp(A dt). Holding a fission power
instead makes the flux -- and so the matrix -- depend on the composition, which
is what the higher-order schemes exist for, and is the only way to compare
their coefficients and the constant-power normalization against OpenMC's.

Approach:
  * The chain, the one-group micro cross sections and the step-0 densities are
    read from the case directory.
  * The material has a volume of 1 cm^3, so number densities in atom/cm^3 are
    the atom counts density.csv holds and the C++ engine marches.
  * The power is the fission power at step 0 under the case's step-1 flux, so
    the trajectory burns at the rate the fixed-flux case did at its start.
  * Each scheme runs OpenMC's IndependentOperator in 'fission-q' mode over the
    case's schedule, in a directory of its own.

Run from the repository root, e.g.:
    python validation/openmc/generate_constant_power_refs.py --case-dir tests/validation/data/vera_pin1a
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import openmc
import openmc.deplete

JOULE_PER_EV = 1.602176634e-19
BARN_CM2 = 1.0e-24

INTEGRATORS = {
    "predictor": openmc.deplete.PredictorIntegrator,
    "cecm": openmc.deplete.CECMIntegrator,
    "celi": openmc.deplete.CELIIntegrator,
    "leqi": openmc.deplete.LEQIIntegrator,
    "cf4": openmc.deplete.CF4Integrator,
}


def read_case(case_dir: Path):
    with open(case_dir / "schedule.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    dts = [float(r["dt_seconds"]) for r in rows]
    flux = float(rows[0]["flux"])

    xs = defaultdict(dict)  # nuclide -> reaction -> barn
    with open(case_dir / "micro_xs.csv", newline="") as f:
        for r in csv.DictReader(f):
            xs[r["nuclide"]][r["reaction"]] = float(r["xs_barn"])

    atoms0 = {}
    with open(case_dir / "density.csv", newline="") as f:
        for r in csv.DictReader(f):
            if int(r["step"]) == 0 and float(r["atoms"]) > 0.0:
                atoms0[r["nuclide"]] = float(r["atoms"])
    return dts, flux, xs, atoms0


def micro_xs_object(xs) -> openmc.deplete.MicroXS:
    nuclides = sorted(xs)
    reactions = sorted({rx for per in xs.values() for rx in per})
    data = np.zeros((len(nuclides), len(reactions), 1))
    for i, nuc in enumerate(nuclides):
        for j, rx in enumerate(reactions):
            data[i, j, 0] = xs[nuc].get(rx, 0.0)
    return openmc.deplete.MicroXS(data, nuclides, reactions)


def fission_power(chain_file: Path, xs, atoms, flux: float) -> float:
    """sum_j n_j sigma_f,j phi Q_j -- the C++ DepletionSystem::powerFor()."""
    chain = openmc.deplete.Chain.from_xml(str(chain_file))
    power = 0.0
    for nuc, n in atoms.items():
        sigma = xs.get(nuc, {}).get("fission", 0.0)
        if sigma == 0.0 or nuc not in chain:
            continue
        q = next((r.Q for r in chain[nuc].reactions if r.type == "fission"), 0.0)
        power += n * sigma * BARN_CM2 * flux * q * JOULE_PER_EV
    return power


def run_scheme(name, chain_file, micro, atoms0, flux, dts, power):
    """One OpenMC march; returns {nuclide: [atoms at step 0..N]}."""
    volume = 1.0  # cm^3: atom/cm^3 == atoms
    densities = {nuc: n / volume * BARN_CM2 for nuc, n in atoms0.items()}  # atom/b-cm
    operator = openmc.deplete.IndependentOperator.from_nuclides(
        volume,
        densities,
        flux=np.array([flux * volume]),
        micro_xs=micro,
        chain_file=str(chain_file),
        nuc_units="atom/b-cm",
        normalization_mode="fission-q",
    )
    integrator = INTEGRATORS[name](operator, dts, power=power, timestep_units="s")
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as work:
        os.chdir(work)
        try:
            integrator.integrate()
            results = openmc.deplete.Results("depletion_results.h5")
            mat = str(operator.materials[0].id)
            out = {}
            for nuc in results[0].index_nuc:
                _, atoms = results.get_atoms(mat, nuc)
                out[nuc] = [float(a) for a in atoms]
        finally:
            os.chdir(cwd)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-dir", default="tests/validation/data/vera_pin1a",
                    help="case written by generate_vera_pin.py")
    ap.add_argument("--schemes", nargs="+", default=list(INTEGRATORS),
                    choices=list(INTEGRATORS))
    args = ap.parse_args()

    case_dir = Path(args.case_dir).resolve()
    chain_file = case_dir / "chain.xml"
    dts, flux, xs, atoms0 = read_case(case_dir)
    micro = micro_xs_object(xs)
    power = fission_power(chain_file, xs, atoms0, flux)

    out_dir = case_dir / "constant_power"
    out_dir.mkdir(exist_ok=True)
    with open(out_dir / "power.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["power_w"])
        w.writerow([repr(power)])

    for name in args.schemes:
        traj = run_scheme(name, chain_file, micro, atoms0, flux, dts, power)
        with open(out_dir / f"{name}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["step", "nuclide", "atoms"])
            for nuc in sorted(traj):
                for step, n in enumerate(traj[nuc]):
                    w.writerow([step, nuc, repr(n)])
        print(f"{name}: wrote {out_dir / (name + '.csv')}")
    print(f"power {power:.6e} W (fission power at step 0 under flux {flux:.3e})")


if __name__ == "__main__":
    main()
