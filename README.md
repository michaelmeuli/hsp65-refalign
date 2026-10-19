# hsp65-refalign

**Question:** does plain reference alignment of the routine Sanger hsp65 reads (`mlsa-kansasii`
`main3_reference_alignment`) identify the *M. kansasii*-complex species of every isolate that has such a
read, once contaminated isolates are excluded?

**Truth:** the whole-genome species (GTDB-Tk) in `output/screening_map_results.csv`. **Reads:** every hsp65
`.ab1` in `data/sanger/seq_kansasii` (made by `mlsa-kansasii/scripts/03_link_sanger_kansasii.sh`), assigned
to its isolate (`PROBENNUMMER`) through all TNR columns, so all reads of an isolate count.

The repo is self-contained: `hsp65_refalign/refalign.py` and `sanger_io.py` are vendored copies of the
mlsa-kansasii modules (same matching logic, no PDF code), and it aligns the reads itself instead of reading
main3's output. Only the shared data (`data/`, `runs/`, `output/screening_map_results.csv`) is needed.

## Run

```bash
sbatch submit_validate.sbatch     # env_mlsa; about 15 min; log lands in the directory you ran sbatch from
python 02_blast_reads.py          # login node is fine; SENDS two read sequences to NCBI BLAST (internet)
python 03_figure.py               # identity and margin per species -> output/hsp65_refalign/figures/ and article/figures/
sbatch submit_korean_check.sbatch # 04: Korean-type hsp65 check on the isolates' own assemblies
```

Environment: `env_mlsa` (`environment.yml`). `04_korean_type_check.py` lays the 441 bp hsp65 amplicon of
the three Korean GTDB kansasii genomes over every isolate's assembly and writes `korean_type_check.tsv` and
`KOREAN_TYPE_SUMMARY.txt` to `output/hsp65_refalign/`.

`01` writes `output/hsp65_refalign/`: `reads.tsv` (one row per read), `isolates.tsv` (one row per isolate with
its verdict), `SUMMARY.txt`, `nonmatching_reads.fa`. `02` writes `blast/` there (raw XML, `blast_hits.tsv`,
`BLAST_SUMMARY.txt`). Paths are in `hsp65_refalign/__init__.py` (root overridable with `KANSASII_ROOT`).

## Verdicts per isolate

An isolate is **correct** if it has at least one read with status `ok` (>= 98% identity, >= 2 fewer
differences than the next species) and all `ok` reads give the genome species. It is **WRONG** if an `ok` read
disagrees. Isolates with no `ok` read are **no call**. To say why, each non-`ok` read is also aligned to the
isolate's own assembly: below 90% identity, the read does not come from that genome. Traces that are empty or
< 100 bp after trimming count as unusable.

## Result (2026-10-16)

See `SUMMARY.txt` for the current numbers. 136 complex isolates have an hsp65 read; 3 are contaminated
(Mkan329-035, -106, -113) and excluded, leaving 133:

- **109 get a call and all 109 are correct** (0 wrong).
- **24 get no call, so "every isolate with an hsp65 read is identified" does not hold:**
  - **13: the read does not match the isolate's own genome** (~80% to every reference and ~75-80% to the
    isolate's assembly) although the trace is clean. These reads fall into two sequence groups
    (a 16-read group over 9 isolates, 2015-2026, 94-100% identical; a pair from Mkan329-007/-178). BLAST of
    one read per group (`02_blast_reads.py`) gives ***Rothia mucilaginosa*** (98.9% and 97.6% over the whole
    read; an oral commensal, so the TB-11/TB-12 primers co-amplify it from the culture or sample). These isolates'
    genomes are clean (CheckM <= 3.3%) and their 16S/other reads match the complex species, so the hsp65 read
    is a contaminant read, not a wrong isolate. Only the two representatives were BLASTed; the other reads are
    assigned to a group by their identity to it.
  - **10: all reads unusable** (empty or very short traces).
  - **1: poor reads** (Mkan329-102, persicum: four reads at 93-94%).
- So the claim holds in the form: **no isolate with a usable hsp65 read of its own culture is misidentified**.
  The routine fix for the 24 isolates without a call is to repeat the PCR (and for the Rothia ones,
  use a mycobacteria-specific primer pair or a fresh subculture).

The earlier version of this repo read main3's output, which was older than the linked reads and silently skips
unusable traces; that gave 125 analysed isolates and 16 no-calls.
