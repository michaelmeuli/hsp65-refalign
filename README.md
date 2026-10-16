# hsp65-refalign

**Question:** does plain reference alignment of the routine Sanger hsp65 reads (`mlsa-kansasii`
`main3_reference_alignment`) identify the *M. kansasii*-complex species of every isolate that has such a
read, once contaminated isolates are excluded?

**Truth:** the whole-genome species (GTDB-Tk) in `output/screening_map_results.csv`. **Reads:** every hsp65
`.ab1` in `data/sanger/seq_kansasii` (made by `mlsa-kansasii/scripts/03_link_sanger_kansasii.sh`), assigned
to its isolate (`PROBENNUMMER`) through all TNR columns, so all reads of an isolate count.

## Run

```bash
# prerequisite: mlsa-kansasii main3 has been run (output/mlsa/main3_reference_alignment/reference_alignment.tsv)
sbatch submit_validate.sbatch      # env_mlsa; about 3 min; log lands in the directory you ran sbatch from
```

Writes `output/hsp65_refalign_validation/`: `reads.tsv` (one row per read), `isolates.tsv` (one row per
isolate with its verdict), `SUMMARY.txt`. It imports `mlsa` from the sibling `mlsa-kansasii` repo.

## Verdicts per isolate

An isolate is **correct** if it has at least one read with main3 status `ok` (>= 98% identity, >= 2 fewer
differences than the next species) and all `ok` reads give the genome species. It is **WRONG** if an `ok` read
disagrees. Isolates with no `ok` read are **no call**. To say why, each non-`ok` read is also aligned to the
isolate's own assembly: below 90% identity, the read does not come from that genome.

## Result (2026-10-16)

See `SUMMARY.txt` in the output directory for the current numbers. Of 127 complex isolates with an hsp65
read, 2 are contaminated (Mkan329-035, -106; excluded), leaving 125:

- **109 get a call and all 109 are correct** (0 wrong; kansasii 91, persicum 11, pseudokansasii 4,
  attenuatum, innocens, ostraviense 1 each).
- **16 get no call, so the claim "every isolate with an hsp65 read is identified" does not hold.** In 13 the reads
  are clean traces that match neither the references (~80%) nor the isolate's own assembly (~75-80%): the
  Sanger sample is not the sequenced isolate (different culture, other organism or mix-up). 2 reads are
  unusable (too short / no alignment) and 1 isolate (Mkan329-102, persicum) has four 93-94% reads that
  never reach `ok`.
- The claim does hold in the form: **no isolate whose hsp65 read matches its own genome is misidentified**.

Mkan329-035 (kansasii + ostraviense mixture) is the only `ok` call that disagrees with the genome species;
it is the excluded contaminated case (`mlsa-kansasii/Mkan329-035_contamination.md`).
