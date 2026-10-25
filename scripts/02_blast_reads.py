#!/usr/bin/env python
"""Identify what the hsp65 reads that do not match their own genome are, by BLASTing one
representative read per sequence group against NCBI core_nt (blastn, not megablast: the reads are ~80%
identical to every kansasii-complex genome).

Needs nonmatching_reads.fa from 01_validate_hsp65_refalign.py and internet access. SENDS THE READ
SEQUENCES TO NCBI. Each query is the longest non-ok read of the named isolate.

Writes to output/hsp65_refalign/blast/: <group>.query.fa, <group>.xml (raw NCBI result) and, for all
groups, blast_hits.tsv (top hits, one row per hit and group) and BLAST_SUMMARY.txt.
Usage: python 02_blast_reads.py            # uses the groups below
       python 02_blast_reads.py --hits 20
"""
from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the hsp65_refalign package

from hsp65_refalign import OUTPUT  # noqa: E402

# group name -> isolate whose longest non-matching read represents the group. Groups come from an
# all-vs-all comparison of the reads: 16 reads of 9 isolates are 94-100% identical (group A),
# Mkan329-007 and Mkan329-178 share another sequence at 96.7% (group B).
GROUPS = {"A_recurring": "Mkan329-164", "B_pair": "Mkan329-178"}


def representative(fasta: str, probe: str) -> tuple[str, str]:
    recs = [r for r in SeqIO.parse(fasta, "fasta") if r.id.startswith(probe + "|")]
    if not recs:
        raise SystemExit(f"no non-matching read of {probe} in {fasta}")
    best = max(recs, key=lambda r: len(r.seq))
    return best.id, str(best.seq)


def parse_hits(xml_path: str, group: str, query_len: int, n: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for hit in ET.parse(xml_path).getroot().iter("Hit"):
        hsp = hit.find("Hit_hsps/Hsp")
        if hsp is None:
            continue
        aln = int(hsp.findtext("Hsp_align-len", "0"))
        rows.append({"group": group, "accession": hit.findtext("Hit_accession"), "title": hit.findtext("Hit_def"),
                     "identity_pct": round(100 * int(hsp.findtext("Hsp_identity", "0")) / max(aln, 1), 2),
                     "aligned_bp": aln, "query_coverage_pct": round(100 * aln / query_len, 1),
                     "bitscore": float(hsp.findtext("Hsp_bit-score", "0")), "evalue": hsp.findtext("Hsp_evalue")})
        if len(rows) == n:
            break
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", type=int, default=10, help="top hits to keep per group (default 10)")
    ap.add_argument("--in-dir", type=Path, default=OUTPUT, help="directory with nonmatching_reads.fa (from 01)")
    ap.add_argument("--out", type=Path, default=OUTPUT / "blast")
    args = ap.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    lines = []
    for group, probe in GROUPS.items():
        name, seq = representative(str(args.in_dir / "nonmatching_reads.fa"), probe)
        (out / f"{group}.query.fa").write_text(f">{name}\n{seq}\n", encoding="utf-8")
        print(f"{group}: BLASTing {name} ({len(seq)} bp) at NCBI ...", flush=True)
        handle = NCBIWWW.qblast("blastn", "core_nt", seq, megablast=False, hitlist_size=max(args.hits, 20),
                                expect=1e-10, format_type="XML")
        xml_path = out / f"{group}.xml"
        xml_path.write_text(handle.read(), encoding="utf-8")
        rows = parse_hits(str(xml_path), group, len(seq), args.hits)
        all_rows += rows
        lines.append(f"{group}: {name} ({len(seq)} bp), top {len(rows)} hits")
        lines += [f"  {r['identity_pct']}% over {r['aligned_bp']} bp ({r['query_coverage_pct']}% of read)  {r['accession']}  {r['title']}"
                  for r in rows] or ["  no hit below e-value 1e-10"]
    pd.DataFrame(all_rows).to_csv(out / "blast_hits.tsv", sep="\t", index=False)
    (out / "BLAST_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
