#!/usr/bin/env python
"""Figure: identity to the closest reference (A) and margin to the next species (B) of every hsp65 read,
per genome species. Reads of all isolates (--exclude-contaminated: uncontaminated only). Reads that are `ok` are coloured by species, the
others (failed, contaminant or poor reads) are grey in A. Dashed lines are the calling thresholds
(98% identity, 2 differences). Needs reads.tsv from 01_validate_hsp65_refalign.py.

Writes figures/hsp65_identity_margin.{png,pdf} to OUTPUT and to article/figures/ in this repo (png).
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the hsp65_refalign package

from hsp65_refalign import OUTPUT  # noqa: E402
from hsp65_refalign.kansasii_paths import SPECIES, SPECIES_COLORS  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e4e4e4"


def jitter(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).uniform(-0.28, 0.28, n)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUTPUT, help="hsp65_refalign output dir (reads.tsv in, figures/ out)")
    ap.add_argument("--exclude-contaminated", action="store_true", help="plot reads of uncontaminated isolates only")
    args = ap.parse_args()
    r = pd.read_csv(args.out / "reads.tsv", sep="\t")
    if args.exclude_contaminated:
        r = r[~r["contaminated"]]
    r = r[r["closest_identity"].notna()]
    species = [s for s in SPECIES if (r["wgs_species"] == s).any()]
    ok = r[r["status"] == "ok"]
    other = r[r["status"] != "ok"]

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    for ax in (a, b):
        ax.set_facecolor("white")
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK)
        ax.set_xlim(-0.6, len(species) - 0.4)

    for i, sp in enumerate(species):
        o = other[other["wgs_species"] == sp]
        a.scatter(i + jitter(len(o), i), o["closest_identity"], s=22, facecolors="none", edgecolors="#9a9a9a", lw=1.0, zorder=2)
        k = ok[ok["wgs_species"] == sp]
        a.scatter(i + jitter(len(k), 10 + i), k["closest_identity"], s=22, color=SPECIES_COLORS[sp], edgecolors="white", lw=0.5, zorder=3)
        b.scatter(i + jitter(len(k), 20 + i), k["margin"], s=22, color=SPECIES_COLORS[sp], edgecolors="white", lw=0.5, zorder=3)
    a.axhline(98, color=MUTED, ls="--", lw=1)
    a.text(len(species) - 0.45, 98.3, "98% identity", ha="right", va="bottom", color=MUTED, fontsize=8)
    b.axhline(2, color=MUTED, ls="--", lw=1)
    b.text(len(species) - 0.45, 2.3, "2 differences", ha="right", va="bottom", color=MUTED, fontsize=8)

    def plural(n: int, word: str) -> str:
        return f"{n} {word}" + ("" if n == 1 else "s")

    labels = [f"{s}\n{plural(ok[ok['wgs_species'] == s]['probennummer'].nunique(), 'isolate')}\n"
              f"{plural(int((ok['wgs_species'] == s).sum()), 'read')}" for s in species]
    for ax in (a, b):
        ax.set_xticks(range(len(species)), labels, fontsize=8)
    a.set_ylabel("Identity to closest reference (%)", color=INK)
    b.set_ylabel("Differences to next-closest species", color=INK)
    b.set_ylim(0, None)
    a.set_title("A  Identity", loc="left", fontsize=10, color=INK, fontweight="bold")
    b.set_title("B  Margin (ok reads)", loc="left", fontsize=10, color=INK, fontweight="bold")
    a.scatter([], [], s=22, facecolors="none", edgecolors="#9a9a9a", label=f"not ok ({len(other)} reads)")
    a.scatter([], [], s=22, color=INK, label=f"ok ({len(ok)} reads)")
    a.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(0.0, 0.55))

    out = args.out / "figures"
    out.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"hsp65_identity_margin.{ext}", dpi=300)
    (REPO / "article" / "figures").mkdir(exist_ok=True)
    shutil.copy(out / "hsp65_identity_margin.png", REPO / "article" / "figures" / "hsp65_identity_margin.png")
    print(f"ok reads {len(ok)}, other {len(other)}; wrote {out}")


if __name__ == "__main__":
    main()
