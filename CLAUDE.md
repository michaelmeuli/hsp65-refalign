# CLAUDE.md

See README.md. Easy to get wrong:

- **Depends on mlsa-kansasii:** imports `mlsa` from `../mlsa-kansasii` (path hack in the script) and reads its main3 output
  `output/mlsa/main3_reference_alignment/reference_alignment.tsv`. Rerun main3 first if the reads change.
- **Environment:** `env_mlsa`. Run through `sbatch submit_validate.sbatch`, not the login node.
- **Type check:** `~/data/conda/envs/env_typecheck/bin/mypy` (strict).
- **Outputs** go to `/shares/sander.imm.uzh/MM/kansasii/output/hsp65_refalign_validation/`, not the repo.
- Contamination = non-empty `contamination_flag` in `screening_map_results.csv`; truth species outside `mlsa.SPECIES` are dropped.
