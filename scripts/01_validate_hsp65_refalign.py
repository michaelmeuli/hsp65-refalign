#!/usr/bin/env python
"""Does plain reference alignment of the routine Sanger hsp65 reads identify the M. kansasii-complex
species of every isolate? Truth = whole-genome (GTDB-Tk) species in output/immensekansasii/results/screening_map_results.csv.
Isolates flagged for contamination stay in the primary result (a laboratory cannot know a culture is mixed); SUMMARY.txt
adds the result without them as a sensitivity analysis, and --exclude-contaminated makes that the primary result.

Every hsp65 .ab1 in data/sanger/seq_kansasii is quality-trimmed, assigned to its isolate (PROBENNUMMER,
via all TNR columns) and aligned to one GTDB representative genome per species; the closest species
is the call (`ok`: >= 98% identity and >= 2 fewer differences than the next species). Reads without an
`ok` call are also aligned to the isolate's own assembly, to tell a read that does not match the
sequenced genome (identity < 90%) from a merely poor one.

Writes to output/hsp65-refalign/results/: reads.tsv, isolates.tsv, SUMMARY.txt, nonmatching_reads.fa
(trimmed reads of the isolates whose reads do not match their own genome; input of 02_blast_reads.py).
About 3 s per read plus 10 s per non-ok read: run as a job (submit_validate.sbatch).
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the hsp65_refalign package

from hsp65_refalign import OUTPUT  # noqa: E402
from hsp65_refalign.kansasii_paths import (ASSEMBLIES, GTDB_REPRESENTATIVES, REFERENCE_STRAINS, SANGER_LINK_DIR,  # noqa: E402
                                           SCREENING_MAP_RESULTS, SPECIES, TNR_COLUMNS_RESULTS)
from hsp65_refalign.refalign import Reference, align_to_reference, identify_read, load_references
from hsp65_refalign.sanger_io import classify_locus, load_trimmed_ab1_mixed  # noqa: E402

# a good read of the isolate's own genome aligns at >= 98%; below this it is not from that genome
OWN_GENOME_MIN_IDENTITY = 90.0
NO_MATCH = "no call: read does not match own genome"


def own_assembly_identity(assemblies: Path, probe: str, seq: str) -> float | None:
    """Identity (%) of the read to the isolate's own assembly, None if no alignment or no assembly."""
    fasta = assemblies / probe / "1_unicycler" / f"{probe}.fasta"
    if not fasta.exists():
        return None
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
        return NO_MATCH
    if own.isna().all():
        return "no call: read unusable"
    return "no call: poor read"


def summarise(title: str, d: pd.DataFrame) -> list[str]:
    called = d[d["n_ok"] > 0]
    lines = [f"{title}: {len(d)} isolates"]
    lines += [f"  {v}: {n}" for v, n in d["verdict"].value_counts().items()]
    lines += [f"  isolates with a call: {len(called)}, correct: {(called['verdict'] == 'correct').sum()}",
              "  by WGS species (correct / analysed): "
              + ", ".join(f"{sp} {(g['verdict'] == 'correct').sum()}/{len(g)}" for sp, g in d.groupby("wgs_species"))]
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtdb-representatives", type=Path, default=GTDB_REPRESENTATIVES, help="one genome per species (references)")
    ap.add_argument("--sanger-dir", type=Path, default=SANGER_LINK_DIR, help="directory of linked .ab1 reads")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--assemblies", type=Path, default=ASSEMBLIES, help="<isolate>/1_unicycler/<isolate>.fasta")
    ap.add_argument("--out", type=Path, default=OUTPUT)
    ap.add_argument("--exclude-contaminated", action="store_true",
                    help="make the result without contamination-flagged isolates the primary one (default: keep them)")
    args = ap.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(args.screening_map_results, dtype=str).fillna("")
    tnr2probe = {r[c]: r["PROBENNUMMER"] for _, r in s.iterrows() for c in TNR_COLUMNS_RESULTS if r.get(c)}
    truth = s.set_index("PROBENNUMMER")
    refs = load_references(args.gtdb_representatives, SPECIES)
    print(f"{len(refs)} references", flush=True)

    rows: list[dict[str, object]] = []
    seqs: dict[str, str] = {}
    for f in sorted(args.sanger_dir.glob("*.ab1")):
        if classify_locus(f.name) != "hsp65":
            continue
        m = re.search(r"(\d{10})", f.name)
        probe = tnr2probe.get(m.group(1), "") if m else ""
        if not probe or probe in REFERENCE_STRAINS or truth.at[probe, "species"] not in SPECIES:
            continue
        loaded = load_trimmed_ab1_mixed(f)
        row: dict[str, object] = {"probennummer": probe, "wgs_species": truth.at[probe, "species"],
                                  "contaminated": truth.at[probe, "contamination_flag"] != "", "read": f.stem}
        if loaded is None:
            row["status"] = "unusable"
            rows.append(row)
            continue
        seq = loaded[0].upper()
        res = identify_read(refs, f.stem, "hsp65", seq)
        row.update(read_bp=len(seq), mixed_fraction=None if loaded[2] is None else round(loaded[2], 4), status=res.status,
                   closest_species=res.closest_species)
        if res.best is not None:
            row.update(closest_identity=round(100 * res.identity(res.best), 2), closest_diffs=res.diffs[res.best])
            if res.runner_up is not None:
                row["margin"] = res.diffs[res.runner_up] - res.diffs[res.best]
        if res.status != "ok":
            row["own_assembly_identity"] = own_assembly_identity(args.assemblies, probe, seq)
            seqs[f.stem] = seq
        print(f"{f.stem}: {row['status']} {row.get('closest_species')} {row.get('closest_identity')}", flush=True)
        rows.append(row)
    ra = pd.DataFrame(rows)
    ra.sort_values(["probennummer", "read"]).to_csv(out / "reads.tsv", sep="\t", index=False)

    iso_rows = []
    for probe, d in ra.groupby("probennummer"):
        species = str(d["wgs_species"].iloc[0])
        ok = d[d["status"] == "ok"]
        iso_rows.append({"probennummer": probe, "wgs_species": species, "contaminated": bool(d["contaminated"].iloc[0]),
                         "n_hsp65_reads": len(d), "n_ok": len(ok), "called_species": ",".join(sorted(set(ok["closest_species"]))),
                         "best_identity": d["closest_identity"].max(),
                         "best_own_assembly_identity": d["own_assembly_identity"].max(), "verdict": verdict(d, species)})
    iso = pd.DataFrame(iso_rows)
    iso.to_csv(out / "isolates.tsv", sep="\t", index=False)

    nm = set(iso[iso["verdict"] == NO_MATCH]["probennummer"])
    with open(out / "nonmatching_reads.fa", "w", encoding="utf-8") as fh:
        for _, r in ra[ra["probennummer"].isin(nm) & (ra["status"] != "ok")].iterrows():
            if r["read"] in seqs:
                fh.write(f">{r['probennummer']}|{r['read']}\n{seqs[r['read']]}\n")

    flagged = iso[iso["contaminated"]]
    clean = iso[~iso["contaminated"]]
    primary = clean if args.exclude_contaminated else iso
    lines = [f"isolates with a WGS species in the complex and >= 1 Sanger hsp65 read: {len(iso)}",
             f"flagged for contamination: {', '.join(flagged['probennummer'])}"
             + (" (excluded)" if args.exclude_contaminated else " (kept in the primary result)"),
             ""]
    lines += summarise("primary result" + (" (flagged isolates excluded)" if args.exclude_contaminated else " (all isolates)"), primary)
    if not args.exclude_contaminated:
        lines += ["", *summarise("sensitivity analysis (flagged isolates excluded)", clean)]
    lines += ["", "flagged isolates:"]
    lines += [f"  {r.probennummer} ({r.wgs_species}, called: {str(r.called_species) or '-'}): {r.n_ok}/{r.n_hsp65_reads} ok reads -> {r.verdict}"
              for r in flagged.itertuples()]
    nc = primary[primary["n_ok"] == 0]
    if len(nc):
        lines += ["", "isolates without a call:"]
        lines += [f"  {r.probennummer} ({r.wgs_species}): {r.n_hsp65_reads} reads, best identity to references {r.best_identity}%, "
                  f"to own assembly {r.best_own_assembly_identity}% -> {r.verdict}" for r in nc.itertuples()]
    (out / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
