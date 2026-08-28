# Validation report

Regenerates a formal comparison of this engine against OpenMC from the committed
reference data, so the claim "we reproduce OpenMC to *x*" is a build artifact
rather than a number quoted in prose.

```bash
cmake -B build -DCRAM_WITH_CHAIN_XML=ON
cmake --build build --target validation-report
# -> build/validation/vera_report.json
# -> build/validation/vera_validation_report.html
```

CI runs the same target on every push (the `validation-report` job) and publishes
both files as an artifact, with the headline numbers written to the run summary.

## Two halves, on purpose

| Step | File | Needs |
|------|------|-------|
| Measure | [vera_report.cpp](vera_report.cpp) | the C++ build, `CRAM_WITH_CHAIN_XML` |
| Render | [render_report.py](render_report.py) | any Python 3 interpreter |

The measurement step replays the case through the real engine and writes the
whole comparison — every nuclide, every step — as JSON. It draws no conclusions
and applies no corrections. The rendering step turns that JSON into a
self-contained HTML page with inline SVG figures.

Splitting them keeps the numbers reusable (the JSON is the machine-readable
record; plot it however you like) and keeps the build's only hard dependency on
the C++ side. `render_report.py` uses the standard library alone — no matplotlib,
no numpy — so a CI runner needs nothing installed beyond Python, and the figures
are drawn against the report's own palette tokens rather than a library's
defaults. If no Python 3 interpreter is found at configure time, the target still
writes the JSON and says the HTML was skipped.

## What the report contains

1. **Purpose and scope** — why the comparison is transport-free, and what that
   does and does not establish.
2. **Setup** — chain size, cross-section provenance, initial condition, the full
   irradiation schedule, and the chain reader's own diagnostics.
3. **Numerics** — integrator, CRAM order, flux normalization, and the definition
   of the reported error.
4. **Results** — summary tiles, the benchmark subset quoted separately from the
   full tracked inventory, and five figures:
   - end-of-life error for every nuclide, ranked
   - error against end-of-life inventory (is it only trace nuclides?)
   - error development across the schedule
   - benchmark-nuclide trajectories against the reference
   - the largest deviations, as trajectories
5. **Full results** — every compared nuclide in a table.

The benchmark subset and the full inventory are always reported side by side.
A curated subset that passes says nothing about the nuclides outside it, and the
figures are deliberately drawn over everything for that reason.

## Options

`vera_report` takes `--data DIR`, `--case NAME`, `--out FILE`, `--tolerance X`
and `--fail-over-tolerance N`. The tolerance is a *reporting* threshold: it sets
where the figures draw their reference line and what the summary counts, but the
tool exits 0 regardless. `--fail-over-tolerance N` turns it into a gate, exiting
1 when more than `N` nuclides exceed it, which is how this would become a CI
check. The `validation-report` job reports without gating today; `vera_pin1a`
currently puts 0 of 181 nuclides over 1e-3 (worst 3.8e-13), so the gate has
nothing to catch until a case is added that does.

Set `SOURCE_DATE_EPOCH` to pin the timestamp and get a byte-identical
regeneration.

## Adding a case

Generate a new case directory with
[../openmc/generate_vera_pin.py](../openmc/generate_vera_pin.py), drop it beside
`vera_pin1a` under `tests/validation/data/`, and pass `--case <name>`. The loader
shared by the report and the regression test is
[../vera_case.hpp](../vera_case.hpp); nothing else needs to change.
