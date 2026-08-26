#!/usr/bin/env python3
"""Render a vera_report JSON document into a self-contained HTML validation report.

Deliberately dependency-free (standard library only). The figures are emitted as
inline SVG rather than through a plotting library so that regenerating the report
in CI needs nothing beyond a Python 3 interpreter, and so the chart tokens match
the report's own palette exactly.

    python render_report.py vera_report.json -o vera_validation_report.html
"""

from __future__ import annotations

import argparse
import html
import json
import math
import sys
from pathlib import Path

# --- design tokens ----------------------------------------------------------
# Light-mode instance of the validated reference palette. The report commits to
# a single (print-oriented) look: the figures are baked SVG, so an automatic
# dark flip would leave axis ink unreadable rather than merely different.

SURFACE = "#fcfcfb"
PAGE = "#f9f9f7"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BORDER = "rgba(11,11,11,0.10)"

S1 = "#2a78d6"  # categorical slot 1 - blue   (this engine)
S2 = "#eb6834"  # categorical slot 2 - orange (OpenMC reference)
S3 = "#1baf7a"  # categorical slot 3 - aqua
CRITICAL = "#d03b3b"
GOOD = "#0ca30c"

# IBM Plex: commissioned for technical communication, which is the register this
# document is in. Serif for headings (a published-benchmark feel), sans for
# running text, mono wherever digits line up. Loaded over the network but never
# depended on: each stack falls back to a system face, so the offline copy in a
# CI artifact is still correctly set, just in the fallback.
FONT_IMPORT = ("@import url('https://fonts.googleapis.com/css2?"
               "family=IBM+Plex+Mono:wght@400;600&"
               "family=IBM+Plex+Sans:wght@400;500;600&"
               "family=IBM+Plex+Serif:wght@500;600&display=swap');")
FONT_BODY = '"IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif'
FONT_HEAD = '"IBM Plex Serif", Georgia, "Times New Roman", serif'
FONT_MONO = '"IBM Plex Mono", ui-monospace, "Cascadia Code", Consolas, monospace'

# Figures are data plates: mono throughout so every tick and label is tabular.
FONT = FONT_MONO

SUPER = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def esc(s) -> str:
    return html.escape(str(s), quote=True)


def fmt_sci(x: float, digits: int = 3) -> str:
    if x == 0:
        return "0"
    if not math.isfinite(x):
        return "n/a"
    return f"{x:.{digits}e}"


# --- a small SVG plotting core ----------------------------------------------


class Scale:
    """Maps a data range onto a pixel range, linearly or in log10."""

    def __init__(self, lo, hi, px0, px1, log=False, invert=False):
        self.log = log
        if log:
            lo = max(lo, 1e-300)
            hi = max(hi, lo * 10)
            self.lo, self.hi = math.log10(lo), math.log10(hi)
        else:
            if hi == lo:
                hi = lo + 1.0
            self.lo, self.hi = lo, hi
        self.px0, self.px1 = px0, px1
        self.invert = invert

    def __call__(self, v):
        if self.log:
            v = math.log10(max(v, 1e-300))
        t = (v - self.lo) / (self.hi - self.lo)
        if self.invert:
            t = 1.0 - t
        return self.px0 + t * (self.px1 - self.px0)

    def ticks(self, target=6):
        if self.log:
            lo, hi = math.floor(self.lo), math.ceil(self.hi)
            step = max(1, int(math.ceil((hi - lo) / target)))
            out = []
            e = lo
            while e <= hi:
                out.append((10.0**e, "10" + str(int(e)).translate(SUPER)))
                e += step
            return out
        span = self.hi - self.lo
        raw = span / max(target, 1)
        mag = 10.0 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
        for m in (1, 2, 2.5, 5, 10):
            if raw <= m * mag:
                step = m * mag
                break
        else:
            step = 10 * mag
        out = []
        v = math.ceil(self.lo / step) * step
        while v <= self.hi + 1e-12:
            label = f"{v:g}"
            out.append((v, label))
            v += step
        return out


class Fig:
    """Accumulates SVG for one figure with a single plot frame."""

    def __init__(self, width, height, pad=(52, 16, 42, 64)):
        # pad = (top, right, bottom, left)
        self.w, self.h = width, height
        self.pt, self.pr, self.pb, self.pl = pad
        self.x0, self.x1 = self.pl, width - self.pr
        self.y0, self.y1 = self.pt, height - self.pb
        self.parts: list[str] = []

    # -- primitives
    def add(self, s):
        self.parts.append(s)

    def text(self, x, y, s, fill=INK2, size=11, anchor="start", weight="400"):
        self.add(
            f'<text x="{x:.1f}" y="{y:.1f}" fill="{fill}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}">{esc(s)}</text>'
        )

    def line(self, x1, y1, x2, y2, stroke, width=1):
        self.add(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"/>'
        )

    def polyline(self, pts, stroke, width=2, opacity=1.0):
        if len(pts) < 2:
            return
        d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        self.add(
            f'<polyline points="{d}" fill="none" stroke="{stroke}" stroke-width="{width}" '
            f'stroke-linejoin="round" stroke-linecap="round" opacity="{opacity}"/>'
        )

    def dot(self, x, y, fill, r=4.0, title=None, ring=SURFACE):
        t = f"<title>{esc(title)}</title>" if title else ""
        self.add(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" '
            f'stroke="{ring}" stroke-width="2">{t}</circle>'
        )

    # -- frame
    def frame(self, sx: Scale, sy: Scale, xlabel, ylabel, title=None, subtitle=None,
              xticks=None, yticks=None):
        if title:
            self.text(self.pl - 40, 20, title, fill=INK, size=14, weight="600")
        if subtitle:
            self.text(self.pl - 40, 36, subtitle, fill=MUTED, size=11)
        for v, lab in (yticks if yticks is not None else sy.ticks()):
            y = sy(v)
            if not (self.y0 - 0.5 <= y <= self.y1 + 0.5):
                continue
            self.line(self.x0, y, self.x1, y, GRID, 1)
            self.text(self.x0 - 8, y + 3.5, lab, fill=MUTED, size=10, anchor="end")
        for v, lab in (xticks if xticks is not None else sx.ticks()):
            x = sx(v)
            if not (self.x0 - 0.5 <= x <= self.x1 + 0.5):
                continue
            self.text(x, self.y1 + 16, lab, fill=MUTED, size=10, anchor="middle")
        self.line(self.x0, self.y1, self.x1, self.y1, AXIS, 1)
        self.text((self.x0 + self.x1) / 2, self.h - 6, xlabel, fill=INK2, size=11, anchor="middle")
        self.add(
            f'<text transform="translate(14,{(self.y0 + self.y1) / 2:.1f}) rotate(-90)" '
            f'fill="{INK2}" font-size="11" text-anchor="middle">{esc(ylabel)}</text>'
        )

    def legend(self, entries, x=None, y=None):
        """entries: list of (label, color, kind) where kind is 'line' or 'dot'."""
        x = self.x1 if x is None else x
        y = self.pt - 14 if y is None else y
        cursor = x
        for label, color, kind in reversed(entries):
            wpx = 7 * len(label) + 26
            cursor -= wpx
            cy = y
            if kind == "line":
                self.line(cursor, cy, cursor + 14, cy, color, 2.5)
            else:
                self.dot(cursor + 7, cy, color, r=4)
            self.text(cursor + 20, cy + 4, label, fill=INK2, size=11)

    def svg(self, label):
        return (
            f'<svg viewBox="0 0 {self.w} {self.h}" width="100%" role="img" '
            f'aria-label="{esc(label)}" font-family=\'{FONT}\'>'
            f'<rect width="{self.w}" height="{self.h}" fill="{SURFACE}"/>'
            + "".join(self.parts)
            + "</svg>"
        )


# --- figures ----------------------------------------------------------------


def place_labels(f: Fig, items, min_gap=12.0, dx=8.0, size=10):
    """Draw point labels, nudging them apart vertically so they cannot overlap.

    The worst nuclides in a ranked plot sit at almost the same x, so labelling
    them naively stacks several strings within a few pixels. Walk them in y order
    and push each one below the last if it is too close, drawing a short leader
    when a label has been moved off its point.
    """
    items = sorted(items, key=lambda it: it[1])
    last_y = -1e9
    for x, y, text in items:
        ly = y if y - last_y >= min_gap else last_y + min_gap
        last_y = ly
        if abs(ly - y) > 1.5:  # leader back to the mark it belongs to
            f.line(x + dx - 2, ly - 3.5, x + 2.5, y, AXIS, 1)
        f.text(x + dx, ly, text, fill=INK, size=size, weight="600")


def fig_error_rank(d) -> str:
    rows = d["nuclides"]
    tol = d["tolerance"]
    floor = 1e-14
    f = Fig(880, 380)
    errs = [max(r["final_relerr"], floor) for r in rows]
    if not errs:
        return ""
    sx = Scale(0, max(len(rows) - 1, 1), f.x0, f.x1)
    sy = Scale(min(errs), max(max(errs), tol * 10), f.y1, f.y0, log=True)
    f.frame(sx, sy, "nuclides, ranked by end-of-life relative error",
            "relative error vs OpenMC",
            title="End-of-life agreement, every tracked nuclide",
            subtitle=f"{len(rows)} nuclides with a positive reference density; "
                     f"log scale, tolerance {fmt_sci(tol, 0)}")
    # tolerance reference line (status colour, always directly labelled)
    ty = sy(tol)
    f.line(f.x0, ty, f.x1, ty, CRITICAL, 1.5)
    f.text(f.x1 - 4, ty - 7, f"tolerance {fmt_sci(tol, 0)}", fill=CRITICAL, size=10,
           anchor="end", weight="600")
    for i, r in enumerate(rows):
        e = max(r["final_relerr"], floor)
        colour = S1 if r["curated"] else MUTED
        f.dot(sx(i), sy(e), colour, r=3.5,
              title=f'{r["name"]}: relative error {fmt_sci(r["final_relerr"])}')
    # direct labels on the worst few, so the outliers are named on the page
    place_labels(f, [(sx(i), sy(max(r["final_relerr"], floor)) + 4, r["name"])
                     for i, r in enumerate(rows[:5])])
    f.legend([("benchmark nuclide", S1, "dot"), ("other tracked nuclide", MUTED, "dot")])
    return f.svg("End-of-life relative error for every tracked nuclide, ranked")


def fig_error_vs_abundance(d) -> str:
    rows = d["nuclides"]
    tol = d["tolerance"]
    last = d["schedule"]["steps"]
    floor = 1e-14
    pts = [(r["ref"][last], max(r["final_relerr"], floor), r) for r in rows if r["ref"][last] > 0]
    if not pts:
        return ""
    f = Fig(880, 380)
    sx = Scale(min(p[0] for p in pts), max(p[0] for p in pts), f.x0, f.x1, log=True)
    sy = Scale(min(p[1] for p in pts), max(max(p[1] for p in pts), tol * 10), f.y1, f.y0, log=True)
    f.frame(sx, sy, "OpenMC end-of-life inventory (atoms)", "relative error vs OpenMC",
            title="Does the disagreement sit in the trace nuclides?",
            subtitle="Each point is one nuclide at end of life")
    ty = sy(tol)
    f.line(f.x0, ty, f.x1, ty, CRITICAL, 1.5)
    f.text(f.x1 - 4, ty - 7, f"tolerance {fmt_sci(tol, 0)}", fill=CRITICAL, size=10,
           anchor="end", weight="600")
    for x, y, r in pts:
        f.dot(sx(x), sy(y), S1 if r["curated"] else MUTED, r=3.5,
              title=f'{r["name"]}: {fmt_sci(x)} atoms, relative error {fmt_sci(r["final_relerr"])}')
    place_labels(f, [(sx(x), sy(y) + 4, r["name"])
                     for x, y, r in sorted(pts, key=lambda p: -p[1])[:5]])
    f.legend([("benchmark nuclide", S1, "dot"), ("other tracked nuclide", MUTED, "dot")])
    return f.svg("Relative error against end-of-life inventory, log-log")


def fig_error_vs_time(d) -> str:
    rows = d["nuclides"]
    days = d["schedule"]["cumulative_days"]
    tol = d["tolerance"]
    n = len(days)
    worst, median, worst_bench = [], [], []
    for k in range(n):
        vals = [r["relerr"][k] for r in rows if r["relerr"][k] is not None]
        bench = [r["relerr"][k] for r in rows if r["curated"] and r["relerr"][k] is not None]
        if not vals:
            worst.append(None); median.append(None); worst_bench.append(None); continue
        vals.sort()
        worst.append(vals[-1])
        median.append(vals[len(vals) // 2])
        worst_bench.append(max(bench) if bench else None)
    floor = 1e-14
    finite = [v for v in worst + median + worst_bench if v]
    if not finite:
        return ""
    f = Fig(880, 360)
    sx = Scale(0, max(days), f.x0, f.x1)
    sy = Scale(max(min(finite), floor), max(max(finite), tol * 10), f.y1, f.y0, log=True)
    f.frame(sx, sy, "cumulative burn time (days)", "relative error vs OpenMC",
            title="How the disagreement develops over the schedule",
            subtitle="Across all tracked nuclides at each reported step")

    ty = sy(tol)
    f.line(f.x0, ty, f.x1, ty, CRITICAL, 1.5)
    f.text(f.x1 - 4, ty - 7, f"tolerance {fmt_sci(tol, 0)}", fill=CRITICAL, size=10,
           anchor="end", weight="600")

    for series, colour, label in ((worst, S2, "worst nuclide"),
                                  (worst_bench, S3, "worst benchmark nuclide"),
                                  (median, S1, "median nuclide")):
        pts = [(sx(days[k]), sy(max(series[k], floor))) for k in range(n) if series[k]]
        f.polyline(pts, colour, 2)
        for k in range(n):
            if series[k]:
                f.dot(sx(days[k]), sy(max(series[k], floor)), colour, r=3.5,
                      title=f"{label} at day {days[k]:.1f}: {fmt_sci(series[k])}")
    f.legend([("worst nuclide", S2, "line"), ("worst benchmark nuclide", S3, "line"),
              ("median nuclide", S1, "line")])
    return f.svg("Relative error against burn time")


def fig_trajectories(d, names, title, subtitle, cols=3) -> str:
    by_name = {r["name"]: r for r in d["nuclides"]}
    picks = [by_name[n] for n in names if n in by_name]
    if not picks:
        return ""
    days = d["schedule"]["cumulative_days"]
    rows_n = (len(picks) + cols - 1) // cols
    pw, ph = 280, 190
    W, H = cols * pw + 20, rows_n * ph + 66
    f = Fig(W, H, pad=(56, 10, 10, 10))
    f.text(10, 20, title, fill=INK, size=14, weight="600")
    f.text(10, 36, subtitle, fill=MUTED, size=11)
    f.legend([("this engine", S1, "line"), ("OpenMC reference", S2, "dot")], x=W - 10, y=22)

    for idx, r in enumerate(picks):
        cx, cy = idx % cols, idx // cols
        ox, oy = 10 + cx * pw, 56 + cy * ph
        px0, px1 = ox + 52, ox + pw - 14
        py0, py1 = oy + 22, oy + ph - 30
        vals = [v for v in r["cram"] + r["ref"] if v > 0]
        if not vals:
            continue
        sx = Scale(0, max(days), px0, px1)
        sy = Scale(min(vals), max(vals), py1, py0, log=True)
        for v, lab in sy.ticks(target=3):
            y = sy(v)
            if py0 - 0.5 <= y <= py1 + 0.5:
                f.line(px0, y, px1, y, GRID, 1)
                f.text(px0 - 6, y + 3.5, lab, fill=MUTED, size=9, anchor="end")
        for v, lab in sx.ticks(target=3):
            x = sx(v)
            if px0 - 0.5 <= x <= px1 + 0.5:
                f.text(x, py1 + 14, lab, fill=MUTED, size=9, anchor="middle")
        f.line(px0, py1, px1, py1, AXIS, 1)
        f.text(px0, oy + 14, r["name"], fill=INK, size=12, weight="600")
        f.text(px1, oy + 14, f'{fmt_sci(r["final_relerr"], 2)}', fill=MUTED, size=10, anchor="end")

        ref_pts = [(sx(days[k]), sy(r["ref"][k])) for k in range(len(days)) if r["ref"][k] > 0]
        f.polyline(ref_pts, S2, 2)
        cram_pts = [(sx(days[k]), sy(r["cram"][k])) for k in range(len(days)) if r["cram"][k] > 0]
        f.polyline(cram_pts, S1, 2)
        for k in range(len(days)):
            if r["ref"][k] > 0:
                f.dot(sx(days[k]), sy(r["ref"][k]), S2, r=3,
                      title=(f'{r["name"]} day {days[k]:.1f} — OpenMC {fmt_sci(r["ref"][k])}, '
                             f'cram {fmt_sci(r["cram"][k])}'))
        f.text((px0 + px1) / 2, oy + ph - 12, "days", fill=MUTED, size=9, anchor="middle")
    return f.svg(title)


# --- HTML -------------------------------------------------------------------

CSS = f"""{FONT_IMPORT}
:root {{ color-scheme: light; }}
* {{ box-sizing: border-box; }}
body {{ margin:0; background:{PAGE}; color:{INK}; font-family:{FONT_BODY};
       font-size:15px; line-height:1.65; -webkit-font-smoothing:antialiased; }}
.wrap {{ max-width: 960px; margin: 0 auto; padding: 48px 24px 96px; }}
h1, h2, h3 {{ font-family:{FONT_HEAD}; text-wrap:balance; }}
h1 {{ font-size: 34px; font-weight:600; line-height:1.15; margin:0 0 8px;
      letter-spacing:-0.015em; }}
h2 {{ font-size: 21px; font-weight:600; margin: 46px 0 12px; padding-top:20px;
      border-top:1px solid {BORDER}; }}
h3 {{ font-size: 16px; font-weight:600; margin: 28px 0 8px; color:{INK}; }}
p, li {{ color:{INK2}; max-width: 68ch; }}
code {{ font-family:{FONT_MONO}; font-size:0.88em;
        background:rgba(11,11,11,0.05); padding:1px 5px; border-radius:4px; }}
.sub {{ color:{MUTED}; font-size:13px; margin:0 0 28px; font-family:{FONT_MONO};
        letter-spacing:-0.01em; }}
.tiles {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px;
          margin:24px 0 8px; }}
.tile {{ background:{SURFACE}; border:1px solid {BORDER}; border-radius:10px; padding:14px 16px; }}
.tile .k {{ font-size:10.5px; text-transform:uppercase; letter-spacing:.09em; color:{MUTED};
            font-weight:500; }}
.tile .v {{ font-size:25px; font-weight:600; line-height:1.25; margin-top:3px;
            font-family:{FONT_MONO}; font-variant-numeric: tabular-nums;
            letter-spacing:-0.02em; }}
.tile .n {{ font-size:12px; color:{MUTED}; margin-top:1px; }}
figure {{ margin:22px 0; background:{SURFACE}; border:1px solid {BORDER}; border-radius:10px;
          padding:8px 8px 4px; overflow-x:auto; }}
figcaption {{ font-size:12.5px; color:{MUTED}; padding:2px 10px 10px; }}
table {{ border-collapse:collapse; width:100%; font-size:12.5px; background:{SURFACE};
         font-family:{FONT_MONO}; font-variant-numeric: tabular-nums; }}
.scroll {{ overflow-x:auto; border:1px solid {BORDER}; border-radius:10px; }}
.tallscroll {{ max-height:520px; overflow-y:auto; }}
th, td {{ text-align:left; padding:7px 12px; border-bottom:1px solid {GRID}; white-space:nowrap; }}
th {{ position:sticky; top:0; background:{SURFACE}; color:{INK}; font-weight:600; font-size:12px;
      text-transform:uppercase; letter-spacing:.05em; z-index:1; }}
td.num, th.num {{ text-align:right; }}
tr:last-child td {{ border-bottom:none; }}
.pill {{ display:inline-block; padding:1px 8px; border-radius:999px; font-size:11.5px;
         font-weight:600; }}
.ok {{ color:{GOOD}; }} .bad {{ color:{CRITICAL}; }}
.badge {{ display:inline-flex; align-items:center; gap:7px; font-weight:600; font-size:13px;
          border:1px solid {BORDER}; border-radius:999px; padding:5px 13px; background:{SURFACE}; }}
.note {{ background:{SURFACE}; border:1px solid {BORDER}; border-left:3px solid {S1};
         border-radius:8px; padding:12px 16px; margin:16px 0; font-size:14px; }}
footer {{ margin-top:48px; padding-top:18px; border-top:1px solid {BORDER};
          color:{MUTED}; font-size:12.5px; }}
"""


def tile(k, v, n="", cls=""):
    return (f'<div class="tile"><div class="k">{esc(k)}</div>'
            f'<div class="v {cls}">{v}</div><div class="n">{esc(n)}</div></div>')


def render(d, fragment: bool = False) -> str:
    s = d["summary"]
    sched = d["schedule"]
    chain = d["chain"]
    tol = d["tolerance"]
    bench = s["benchmark_subset"]
    over = s["over_tolerance"]
    passing = over == 0

    verdict = (f'<span class="badge"><span class="ok">&#10003;</span> '
               f'All {s["compared"]} nuclides within {fmt_sci(tol, 0)}</span>') if passing else (
        f'<span class="badge"><span class="bad">&#9679;</span> '
        f'{over} of {s["compared"]} nuclides exceed {fmt_sci(tol, 0)}</span>')

    # schedule table
    sched_rows = "".join(
        f'<tr><td class="num">{i + 1}</td><td class="num">{sched["dt_seconds"][i]:,.0f}</td>'
        f'<td class="num">{sched["dt_seconds"][i] / 86400.0:.2f}</td>'
        f'<td class="num">{sched["flux"][i]:.4g}</td>'
        f'<td class="num">{sched["cumulative_days"][i + 1]:.2f}</td></tr>'
        for i in range(sched["steps"]))

    # per-nuclide table
    def row(r):
        cls = "bad" if r["final_relerr"] > tol else ""
        mark = "&#9679; " if r["curated"] else ""
        return (f'<tr><td>{mark}{esc(r["name"])}</td>'
                f'<td class="num">{fmt_sci(r["cram"][-1])}</td>'
                f'<td class="num">{fmt_sci(r["ref"][-1])}</td>'
                f'<td class="num {cls}">{fmt_sci(r["final_relerr"])}</td></tr>')

    nuc_rows = "".join(row(r) for r in d["nuclides"])

    worst_names = [r["name"] for r in d["nuclides"][:6]]
    no_yields = chain["fission_without_yield_table"]
    sf_only = chain["decay_only_spontaneous_fission"]

    structural = ""
    if no_yields or sf_only:
        bits = []
        if no_yields:
            bits.append(
                f'<p><strong>{len(no_yields)} fission channel(s) have no usable yield table '
                f'in this chain.</strong> The chain defines a <code>fission</code> reaction for '
                f'these nuclides, but no fission-yield set is reachable for them, so their '
                f'fission products are not determined by the input: '
                f'<code>{esc(", ".join(no_yields))}</code>.</p>')
        if sf_only:
            bits.append(
                f'<p><strong>{len(sf_only)} nuclide(s) decay only by spontaneous fission.</strong> '
                f'Every decay mode the chain gives them is a fission mode, so their entire decay '
                f'is routed through the fission-yield tables: '
                f'<code>{esc(", ".join(sf_only))}</code>.</p>')
        structural = ('<h3>Under-determined topology in this chain</h3>'
                      '<div class="note">' + "".join(bits) +
                      '<p style="margin-bottom:0">These are properties of the input chain, '
                      'recorded because they identify exactly which nuclides the reference '
                      'topology leaves under-constrained. Compare against the outliers above '
                      'before attributing a deviation to the solver.</p></div>')

    diag = chain["diagnostics"]
    diag_rows = "".join(
        f'<tr><td>{esc(k)}</td><td class="num">{v}</td></tr>' for k, v in (
            ("nuclides that did not parse", diag["unparsed_nuclides"]),
            ("decay targets that did not parse", diag["unparsed_decay_targets"]),
            ("yield products that did not parse", diag["unparsed_yield_products"]),
            ("reactions of an unmodeled type", diag["unmodeled_reactions"]),
        ))

    body = f"""<div class="wrap">

<h1>Depletion validation against OpenMC</h1>
<p class="sub">Case <strong>{esc(d["case"]["name"])}</strong> &middot;
cram {esc(d["cram_version"])} &middot; generated {esc(d["generated_utc"])}</p>
<p>{verdict}</p>

<div class="tiles">
{tile("Nuclides compared", s["compared"], "positive reference at end of life")}
{tile("Over tolerance", f'{over}', f'threshold {fmt_sci(tol, 0)}', "ok" if passing else "bad")}
{tile("Worst error", fmt_sci(s["worst"], 3), esc(s["worst_nuclide"]))}
{tile("Median error", fmt_sci(s["median"], 3), "across all compared nuclides")}
</div>

<h2>1. Purpose and scope</h2>
<p>This repository is a transport-free depletion <em>engine</em>. The full VERA depletion
benchmark couples neutron transport to the depletion solver, and transport supplies the
one-group cross sections the burnup matrix is built from &mdash; which this engine does not
compute. The comparison here therefore removes transport from the loop entirely.</p>
<p>OpenMC performed one beginning-of-life transport solve to produce a <strong>fixed</strong>
one-group micro-cross-section set, then depleted at a fixed flux with its first-order
predictor. Given that same chain, the same cross sections and the same flux schedule, this
engine must reproduce OpenMC's number densities to CRAM precision. What is under test is
therefore <strong>burnup-matrix assembly and the matrix exponential</strong>, isolated from
the transport and cross-section error a transport-free engine cannot reproduce.</p>
<p>Agreement here is a statement about this engine versus OpenMC on identical inputs. It is
not a statement about either code versus measurement.</p>

<h2>2. Setup</h2>
<h3>Inputs</h3>
<div class="scroll"><table>
<tr><th>Input</th><th>Value</th></tr>
<tr><td>Case directory</td><td><code>{esc(d["case"]["path"])}</code></td></tr>
<tr><td>Chain</td><td>{chain["nuclides"]} nuclides,
{chain["reaction_channels"]} reaction channels (OpenMC <code>depletion_chain</code> XML)</td></tr>
<tr><td>Cross sections</td><td>fixed one-group micro xs, applied unchanged at every step</td></tr>
<tr><td>Initial condition</td><td>{esc(d["solver"]["initial_condition"])}</td></tr>
<tr><td>Depletion steps</td><td>{sched["steps"]}</td></tr>
<tr><td>Total burn time</td><td>{sched["total_days"]:.2f} days
({sched["total_seconds"]:,.0f} s)</td></tr>
</table></div>

<h3>Chain reader diagnostics</h3>
<p>Entries present in the chain file that the reader could not use. All four should be zero
for a reference case; a non-zero count means the topology under test is not the topology in
the file.</p>
<div class="scroll"><table>
<tr><th>Dropped entries</th><th class="num">Count</th></tr>
{diag_rows}
</table></div>

<h3>Irradiation schedule</h3>
<div class="scroll tallscroll"><table>
<tr><th class="num">Step</th><th class="num">dt (s)</th><th class="num">dt (days)</th>
<th class="num">flux (n/cm&sup2;/s)</th><th class="num">cumulative (days)</th></tr>
{sched_rows}
</table></div>

<h2>3. Numerics</h2>
<div class="scroll"><table>
<tr><th>Setting</th><th>Value</th></tr>
<tr><td>Time integrator</td><td>{esc(d["solver"]["integrator"])} (matches the scheme OpenMC
ran, so integrator error cancels rather than being compared)</td></tr>
<tr><td>Matrix exponential</td><td>CRAM order {d["solver"]["cram_order"]}</td></tr>
<tr><td>Flux normalization</td><td>{esc(d["solver"]["normalization"])}</td></tr>
<tr><td>Reporting tolerance</td><td>{fmt_sci(tol, 0)} on relative error</td></tr>
</table></div>
<p>The error reported throughout is the relative deviation of this engine's inventory from
OpenMC's, per nuclide and per step:</p>
<p style="text-align:center"><code>relerr(i, k) = |n_cram(i, k) &minus; n_openmc(i, k)| /
n_openmc(i, k)</code></p>
<p>A nuclide is compared at a step only where OpenMC reports a strictly positive inventory
for it; a zero reference has no relative error to speak of. Because both codes march the
same integrator over the same fixed matrix, the residual is dominated by matrix assembly
and the exponential, not by time discretization.</p>

<h2>4. Results</h2>
<div class="tiles">
{tile("Benchmark nuclides", f'{bench["compared"]}', "the reported set of interest")}
{tile("&mdash; over tolerance", f'{bench["over_tolerance"]}',
      "within the benchmark subset", "ok" if bench["over_tolerance"] == 0 else "bad")}
{tile("&mdash; worst", fmt_sci(bench["worst"], 3), "benchmark subset only")}
{tile("95th percentile", fmt_sci(s["p95"], 3), "all compared nuclides")}
</div>
<p>The benchmark subset is the actinide and fission-product set the VERA benchmark reports.
It is quoted separately from the full tracked inventory because the two can disagree
sharply: a curated subset that passes says nothing about the nuclides outside it, and the
figures below are drawn over <strong>all</strong> {s["compared"]} compared nuclides for
exactly that reason.</p>

<figure>{fig_error_rank(d)}
<figcaption>Figure 1 &mdash; End-of-life relative error for every nuclide with a positive
reference inventory, ranked worst-first. The five worst are labelled. Benchmark nuclides are
highlighted; hover any point for its value.</figcaption></figure>

<figure>{fig_error_vs_abundance(d)}
<figcaption>Figure 2 &mdash; The same errors against end-of-life inventory. Points to the
left are trace nuclides. This separates &ldquo;small absolute discrepancy in a rare
nuclide&rdquo; from a discrepancy in the bulk inventory.</figcaption></figure>

<figure>{fig_error_vs_time(d)}
<figcaption>Figure 3 &mdash; Error development across the schedule. A flat trace indicates a
constant offset; growth indicates an accumulating one. The worst-nuclide and worst-benchmark
traces are separated so a curated-subset result cannot mask the tail.</figcaption></figure>

<figure>{fig_trajectories(d, list(d["benchmark_nuclides"])[:9],
    "Benchmark nuclide trajectories",
    "This engine against the OpenMC reference; log inventory against burn time")}
<figcaption>Figure 4 &mdash; Trajectories for the benchmark nuclides. The two traces should
overlay; visible separation is the error in Figure 1. End-of-life relative error is printed
at the top right of each panel.</figcaption></figure>

<figure>{fig_trajectories(d, worst_names,
    "Largest deviations",
    "The six nuclides with the greatest end-of-life disagreement")}
<figcaption>Figure 5 &mdash; The six worst-agreeing nuclides. Where a trace runs
systematically above the reference, this engine is under-removing that nuclide.
</figcaption></figure>

{structural}

<h2>5. Full results</h2>
<p>Every nuclide with a positive end-of-life reference inventory, worst first.
&#9679; marks a benchmark nuclide; values over tolerance are shown in red.</p>
<div class="scroll tallscroll"><table>
<tr><th>Nuclide</th><th class="num">this engine (atoms)</th>
<th class="num">OpenMC (atoms)</th><th class="num">relative error</th></tr>
{nuc_rows}
</table></div>

<footer>
<p>Regenerate with <code>cmake --build &lt;build&gt; --target validation-report</code>
(needs <code>-DCRAM_WITH_CHAIN_XML=ON</code>). The measurement step is
<code>validation/report/vera_report.cpp</code>, which writes the JSON this page is rendered
from; the rendering step is <code>validation/report/render_report.py</code> and needs only a
Python 3 interpreter. Reference data is produced offline by
<code>validation/openmc/generate_vera_pin.py</code>.</p>
<p>Report schema {d["schema"]} &middot; generated {esc(d["generated_utc"])} &middot;
cram {esc(d["cram_version"])}</p>
</footer>
</div>"""

    # Named after the case, not captioned: this page sits beside other reports and
    # the case is what tells them apart.
    title = f'{esc(d["case"]["name"])} Depletion Validation'
    if fragment:
        # Body-only: the Artifact host supplies the document skeleton, so emitting
        # our own <html>/<head>/<body> would nest a second document inside it.
        return f"<title>{title}</title>\n<style>{CSS}</style>\n{body}\n"
    return (
        '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{title}</title>\n<style>{CSS}</style></head><body>{body}"
        "</body></html>\n"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("json", type=Path, help="vera_report JSON document")
    ap.add_argument("-o", "--out", type=Path, default=Path("vera_validation_report.html"))
    ap.add_argument("--fragment", action="store_true",
                    help="emit body-only HTML for embedding in a host page")
    args = ap.parse_args(argv)

    if not args.json.exists():
        print(f"render_report: no such file: {args.json}", file=sys.stderr)
        return 2
    with args.json.open(encoding="utf-8") as fh:
        d = json.load(fh)
    if d.get("schema") != 1:
        print(f"render_report: unsupported schema {d.get('schema')!r}", file=sys.stderr)
        return 2

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render(d, fragment=args.fragment), encoding="utf-8")
    s = d["summary"]
    print(f"render_report: {args.out}  "
          f"({s['compared']} nuclides, {s['over_tolerance']} over {d['tolerance']:.0e})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
