#!/usr/bin/env python
"""Does any isolate carry the Korean-type hsp65 (kansasii by genome, persicum-like hsp65)?

The 441 bp hsp65 amplicon of the three Korean GTDB kansasii genomes (GCF_002705785.1, -825.1, -865.1; identical to
each other) is laid directly over every isolate's own Unicycler assembly, so the answer comes from the NGS data and
needs no primers and no Sanger read. For comparison the same is done with the amplicon of the kansasii reference
(ATCC 12478, GCF_000157895.3) and of a persicum reference (GCF_002086675.1).

An isolate is Korean-type if its genome is kansasii and its hsp65 is closer to the Korean amplicon than to the
kansasii reference. Writes korean_type_check.tsv and KOREAN_TYPE_SUMMARY.txt to OUTPUT.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the hsp65_refalign package

from hsp65_refalign import OUTPUT  # noqa: E402
from hsp65_refalign.kansasii_paths import (ASSEMBLIES, ATYPICAL_HSP65_GENOMES, ATCC12478_ACCESSION, HSP65_AMPLICONS,  # noqa: E402
                                           REFERENCE_STRAINS, SCREENING_MAP_RESULTS, SPECIES)
from hsp65_refalign.refalign import Reference, align_to_reference  # noqa: E402

# the first of the three Korean genomes stands for them (near-clonal, ANI 99.9-100); ATCC 12478; M. persicum AFPC-000227
PROBES = {"korean": ATYPICAL_HSP65_GENOMES[0], "kansasii_ref": ATCC12478_ACCESSION, "persicum_ref": "GCF_002086675.1"}


def read_amplicon(amplicons: Path, accession: str) -> str:
    return str(next(SeqIO.parse(amplicons / f"{accession}_hsp65.fasta", "fasta")).seq).upper()


def diffs_to_assembly(fasta: Path, probe: str) -> tuple[int, int] | None:
    """(differences, aligned columns) of the amplicon over the assembly, None if it does not place."""
    hit = align_to_reference(Reference("own", "own", "", fasta), probe)
    if hit is None:
        return None
    cols = [(a, b) for a, b in zip(hit.ref_cols, probe) if a in "ACGT"]
    return sum(a != b for a, b in cols), len(cols)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hsp65-amplicons", type=Path, default=HSP65_AMPLICONS)
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--assemblies", type=Path, default=ASSEMBLIES)
    ap.add_argument("--out", type=Path, default=OUTPUT)
    args = ap.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    amplicons = {k: read_amplicon(args.hsp65_amplicons, a) for k, a in PROBES.items()}
    kor_vs_ref = sum(a != b for a, b in zip(amplicons["korean"], amplicons["kansasii_ref"]))
    print(f"Korean amplicon vs kansasii reference: {kor_vs_ref} differences", flush=True)
    s = pd.read_csv(args.screening_map_results, dtype=str).fillna("")
    rows: list[dict[str, object]] = []
    for _, t in s[s["species"].isin(SPECIES) & ~s["PROBENNUMMER"].isin(REFERENCE_STRAINS)].iterrows():
        pn = t["PROBENNUMMER"]
        fasta = args.assemblies / pn / "1_unicycler" / f"{pn}.fasta"
        if not fasta.exists():
            continue
        row: dict[str, object] = {"probennummer": pn, "wgs_species": t["species"], "contaminated": t["contamination_flag"] != ""}
        for name, seq in amplicons.items():
            r = diffs_to_assembly(fasta, seq)
            row[f"diffs_{name}"], row[f"columns_{name}"] = r if r else (None, 0)
        print(pn, row["diffs_korean"], row["diffs_kansasii_ref"], row["diffs_persicum_ref"], flush=True)
        rows.append(row)
    df = pd.DataFrame(rows)
    df["korean_type"] = df["diffs_korean"] < df["diffs_kansasii_ref"]
    df.to_csv(out / "korean_type_check.tsv", sep="\t", index=False)

    k = df[(df["wgs_species"] == "kansasii") & ~df["contaminated"]]
    lines = [f"Korean amplicon vs kansasii reference amplicon: {kor_vs_ref} differences",
             f"isolates analysed: {len(df)} ({int(df['diffs_korean'].isna().sum())} without placement)",
             f"uncontaminated kansasii isolates: {len(k)}",
             f"  minimum differences to the Korean amplicon: {k['diffs_korean'].min():.0f}",
             f"  maximum differences to the kansasii reference: {k['diffs_kansasii_ref'].max():.0f}",
             f"  closer to the Korean amplicon than to the kansasii reference: {int(k['korean_type'].sum())}",
             f"all isolates closer to the Korean amplicon than to the kansasii reference: {int(df['korean_type'].sum())}"]
    for sp, d in df[~df["contaminated"]].groupby("wgs_species"):
        lines.append(f"  {sp}: n={len(d)}, min diffs to Korean {d['diffs_korean'].min():.0f}, "
                     f"to persicum reference {d['diffs_persicum_ref'].min():.0f}")
    (out / "KOREAN_TYPE_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
