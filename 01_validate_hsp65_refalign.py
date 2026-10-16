#!/usr/bin/env python
"""Does plain reference alignment (mlsa-kansasii main3) of the routine Sanger hsp65 reads identify
the M. kansasii-complex species of every isolate? Truth = whole-genome (GTDB-Tk) species in
output/screening_map_results.csv, contaminated isolates excluded.

Reads main3's output/mlsa/main3_reference_alignment/reference_alignment.tsv (run it first), assigns every
hsp65 read of data/sanger/seq_kansasii to its isolate (PROBENNUMMER, via all TNR columns) and judges
per isolate. Reads without an `ok` call are aligned to the isolate's own assembly to tell a read that
does not match the sequenced genome (low identity, e.g. a different culture or organism) from a merely poor one.

Writes to output/hsp65_refalign_validation/: reads.tsv, isolates.tsv, SUMMARY.txt.
Needs the per-sample assemblies of runs/mkan329/assembly/results; about 10 s per checked read, so run it
as a job (submit_validate.sbatch).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
MLSA_REPO = ROOT / "repos" / "mlsa-kansasii"
sys.path.insert(0, str(MLSA_REPO))

import pandas as pd  # noqa: E402

from mlsa import SANGER_LINK_DIR, SPECIES  # noqa: E402
from mlsa.refalign import Reference, align_to_reference  # noqa: E402
from mlsa.sanger_io import classify_locus, load_trimmed_ab1  # noqa: E402

OUTPUT = ROOT / "output" / "hsp65_refalign_validation"
MAIN3 = ROOT / "output" / "mlsa" / "main3_reference_alignment" / "reference_alignment.tsv"
SCREENING = ROOT / "output" / "screening_map_results.csv"
ASSEMBLIES = ROOT / "runs" / "mkan329" / "assembly" / "results"
TNR_COLUMNS = ["TNR", "TNR_JC", "TNR_MM", "TNR3", "TNR4", "TNR5", "TNR6", "TNR_MLSA"]
# a good read of the isolate's own genome aligns at >= 98%; below this it is not from that genome
OWN_GENOME_MIN_IDENTITY = 90.0


def own_assembly_identity(probe: str, read: str) -> float | None:
    """Identity (%) of the trimmed .ab1 read to the isolate's own assembly, None if no alignment or no assembly."""
    fasta = ASSEMBLIES / probe / "1_unicycler" / f"{probe}.fasta"
    loaded = load_trimmed_ab1(SANGER_LINK_DIR / f"{read}.ab1")
    if not fasta.exists() or loaded is None:
        return None
    seq = loaded[0].upper()
    hit = align_to_reference(Reference("own", "own", "", fasta), seq)
    if hit is None:
        return None
    return round(100 * sum(a == b for a, b in zip(hit.ref_cols, seq)) / len(seq), 2)


def verdict(d: pd.DataFrame, species: str) -> str:
    ok = d[d["status"] == "ok"]
    if len(ok):
        return "correct" if (ok["closest_species"] == species).all() else "WRONG"
    own = d["own_assembly_identity"]
    if own.notna().any() and own.max() < OWN_GENOME_MIN_IDENTITY:
        return "no call: read does not match own genome"
    if own.isna().all():
        return "no call: read unusable"
    return "no call: poor read"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(SCREENING, dtype=str).fillna("")
    tnr2probe = {r[c]: r["PROBENNUMMER"] for _, r in s.iterrows() for c in TNR_COLUMNS if r.get(c)}
    truth = s.set_index("PROBENNUMMER")

    ra = pd.read_csv(MAIN3, sep="\t", dtype={"locus": str}).drop_duplicates("read")
    ra = ra[ra["read"].map(lambda n: classify_locus(n + ".ab1") == "hsp65")].copy()
    ra["probennummer"] = ra["read"].map(lambda n: tnr2probe.get(m.group(1), "") if (m := re.search(r"(\d{10})", n)) else "")
    ra = ra[ra["probennummer"] != ""]
    ra["wgs_species"] = ra["probennummer"].map(truth["species"])
    ra["contaminated"] = ra["probennummer"].map(truth["contamination_flag"]) != ""
    ra = ra[ra["wgs_species"].isin(SPECIES)]
    ra["own_assembly_identity"] = [None if r["status"] == "ok" else own_assembly_identity(r["probennummer"], r["read"])
                                   for _, r in ra.iterrows()]
    cols = ["probennummer", "wgs_species", "contaminated", "read", "read_bp", "mixed_fraction", "status", "closest_species",
            "closest_identity", "closest_diffs", "margin", "own_assembly_identity"]
    ra[cols].sort_values(["probennummer", "read"]).to_csv(OUTPUT / "reads.tsv", sep="\t", index=False)

    rows = []
    for probe, d in ra.groupby("probennummer"):
        species = str(d["wgs_species"].iloc[0])
        ok = d[d["status"] == "ok"]
        rows.append({"probennummer": probe, "wgs_species": species, "contaminated": bool(d["contaminated"].iloc[0]),
                     "n_hsp65_reads": len(d), "n_ok": len(ok),
                     "called_species": ",".join(sorted(set(ok["closest_species"]))),
                     "best_identity": d["closest_identity"].max(), "best_own_assembly_identity": d["own_assembly_identity"].max(),
                     "verdict": verdict(d, species)})
    iso = pd.DataFrame(rows)
    iso.to_csv(OUTPUT / "isolates.tsv", sep="\t", index=False)

    clean = iso[~iso["contaminated"]]
    lines = [f"isolates with a WGS species in the complex and >= 1 Sanger hsp65 read: {len(iso)}",
             f"contaminated (excluded): {', '.join(iso[iso['contaminated']]['probennummer'])}",
             f"analysed: {len(clean)}", ""]
    lines += [f"  {v}: {n}" for v, n in clean["verdict"].value_counts().items()]
    called = clean[clean["n_ok"] > 0]
    lines += ["", f"isolates with a call: {len(called)}, correct: {(called['verdict'] == 'correct').sum()}",
              "by WGS species (correct / analysed): "
              + ", ".join(f"{sp} {(g['verdict'] == 'correct').sum()}/{len(g)}" for sp, g in clean.groupby("wgs_species"))]
    nc = clean[clean["n_ok"] == 0]
    if len(nc):
        lines += ["", "isolates without a call:"]
        lines += [f"  {r.probennummer} ({r.wgs_species}): {r.n_hsp65_reads} reads, best identity to references {r.best_identity}%, "
                  f"to own assembly {r.best_own_assembly_identity}% -> {r.verdict}" for r in nc.itertuples()]
    (OUTPUT / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
