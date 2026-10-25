#!/usr/bin/env python
"""Which GTDB-kansasii assemblies of the public BioProjects carry a persicum-like (Korean-type) hsp65?

The hsp65 amplicon is cut in silico from each sample's own Unicycler assembly (TB-11 / TB-12w, TB-12 as fallback) and
classified like a Sanger read against the 72 per-genome hsp65 amplicons of the GTDB complex genomes, so the three kansasii
genomes with an atypical (persicum-like) hsp65 are separate references. Only samples that GTDB-Tk calls Mycobacterium
kansasii are analysed (--species persicum gives the positive control). Per sample: closest species / accession,
atypical (closest reference is an atypical-hsp65 kansasii genome), persicum_like (closest species is persicum or atypical).
Writes insilico_public.tsv and SUMMARY.txt to OUTPUT_PUBLIC.

Usage: 05_call_public.py [--limit N] [--species S] [--out DIR] [BioProject ...]   (default: the four public projects)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from hsp65_refalign import HSP65_AMPLICONS, OUTPUT_PUBLIC, PUBLIC_PROJECTS, RUNS, SPECIES
from hsp65_refalign.insilico import in_silico_pcr
from hsp65_refalign.refalign import identify_read, load_references

GTDB_KANSASII = "s__Mycobacterium kansasii"


def gtdb_species(project: str) -> dict[str, str]:
    """sample -> GTDB-Tk species of every sample of the project that has a call."""
    out: dict[str, str] = {}
    results = RUNS / f"public_{project}" / "assembly" / "results"
    for f in sorted(results.glob("*/3_quality/GTDB/*_gtdb_summary.tsv")):
        t = pd.read_csv(f, sep="\t", dtype=str)
        for _, r in t.iterrows():
            out[r["user_genome"]] = r["classification"].split(";")[-1]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("projects", nargs="*", default=PUBLIC_PROJECTS)
    ap.add_argument("--limit", type=int, help="only the first N samples per project (test)")
    ap.add_argument("--species", default=GTDB_KANSASII, help="GTDB species to analyse (default kansasii)")
    ap.add_argument("--out", type=Path, default=OUTPUT_PUBLIC)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    refs = load_references(HSP65_AMPLICONS, SPECIES)
    rows: list[dict[str, object]] = []
    for proj in args.projects:
        species = gtdb_species(proj)
        todo = [s for s, sp in species.items() if sp == args.species]
        print(f"{proj}: {len(species)} with GTDB call, {len(todo)} {args.species}", flush=True)
        for sid in todo[:args.limit]:
            f = RUNS / f"public_{proj}" / "assembly" / "results" / sid / "1_unicycler" / f"{sid}.fasta"
            row: dict[str, object] = {"sample": sid, "bioproject": proj, "gtdb_species": species[sid], "amplicon_bp": 0}
            if f.exists():
                amp = in_silico_pcr({r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")})
                if amp:
                    row["amplicon_bp"] = len(amp)
                    res = identify_read(refs, sid, "hsp65", amp)
                    b, s = res.best, res.runner_up
                    row.update(status=res.status, closest_species=res.closest_species, atypical=res.atypical)
                    if b is not None:
                        row.update(closest_accession=res.hits[b].ref.accession, diffs=res.diffs[b],
                                   identity=round(100 * res.identity(b), 2))
                    if b is not None and s is not None:
                        row.update(second_species=res.hits[s].ref.species, margin=res.diffs[s] - res.diffs[b])
            rows.append(row)
            print(row, flush=True)

    df = pd.DataFrame(rows)
    for c in ("status", "closest_species", "atypical", "second_species"):
        if c not in df:
            df[c] = None
    df["persicum_like"] = (df["closest_species"] == "persicum") | (df["atypical"] == True)  # noqa: E712
    df.to_csv(args.out / "insilico_public.tsv", sep="\t", index=False)

    found = df[df["amplicon_bp"] > 0]
    lines = [f"{args.species} samples analysed: {len(df)}; hsp65 amplicon found in {len(found)}, no amplicon in {len(df) - len(found)}", ""]
    lines.append("closest hsp65 reference (all samples with amplicon):")
    for label, mask in (("kansasii (typical)", (found["closest_species"] == "kansasii") & (found["atypical"] != True)),  # noqa: E712
                        ("kansasii, atypical (persicum-like) hsp65 genome", found["atypical"] == True),  # noqa: E712
                        ("persicum", found["closest_species"] == "persicum")):
        lines.append(f"  {label}: {mask.sum()}")
    other = found[~found["closest_species"].isin(["kansasii", "persicum"]) & (found["atypical"] != True)]  # noqa: E712
    lines.append(f"  other species: {len(other)} " + str(other["closest_species"].value_counts().to_dict()))
    pl = found[found["persicum_like"]]
    lines += ["", f"persicum-like hsp65: {len(pl)} of {len(found)}",
              "  by status: " + str(pl["status"].value_counts().to_dict()),
              "  by project: " + str(pl["bioproject"].value_counts().to_dict()),
              f"  kansasii closest but persicum runner-up within 2 diffs (ambiguous): "
              f"{((found['closest_species'] == 'kansasii') & (found['second_species'] == 'persicum') & (found['status'] == 'ambiguous')).sum()}"]
    (args.out / "SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
