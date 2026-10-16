# CLAUDE.md

See README.md. Easy to get wrong:

- **Independent of mlsa-kansasii:** `hsp65_refalign/refalign.py` and `sanger_io.py` are vendored copies (PDF code removed). If the
  matching logic changes there, sync by hand. Do not import `mlsa`.
- **Environment:** `env_mlsa`. Run through `sbatch submit_validate.sbatch`, not the login node.
- **Type check:** `~/data/conda/envs/env_typecheck/bin/mypy` (strict).
- **Outputs** go to `/shares/sander.imm.uzh/MM/kansasii/output/hsp65_refalign/`, not the repo.
- Contamination = non-empty `contamination_flag` in `screening_map_results.csv`; truth species outside `hsp65_refalign.SPECIES` are dropped.
- `02_blast_reads.py` sends sequences to NCBI; only run it when that is intended. Its two groups (`GROUPS`) were chosen by an all-vs-all comparison of the non-matching reads.
