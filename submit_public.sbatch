#!/usr/bin/env bash
# 05: persicum-like hsp65 among GTDB-kansasii assemblies of the public BioProjects -> output/hsp65_refalign/public/
#   sbatch submit_public.sbatch        the real run
#   sbatch submit_public.sbatch persicum  all GTDB-persicum samples (80) -> output/hsp65_refalign/public_persicum/ (hsp65 of the
#                                         samples whose gyrA is kansasii-like; read by kansasii-docs/bin/make_atypical_table.py)
#   sbatch submit_public.sbatch ctl    positive control: first 3 GTDB-persicum samples per project -> output/hsp65_refalign/public_ctl/
#SBATCH --job-name=hsp65-public
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/hsp65_refalign/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/hsp65-refalign
OUTPUT=$OUTPUT_ROOT/hsp65_refalign
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

cd "$REPO"
if [ "${1:-}" = "persicum" ]; then
  OUT="$OUTPUT/public_persicum"
  python scripts/05_call_public.py --species "s__Mycobacterium persicum" --out "$OUT"
elif [ "${1:-}" = "ctl" ]; then
  OUT="$OUTPUT/public_ctl"
  python scripts/05_call_public.py --limit 3 --species "s__Mycobacterium persicum" --out "$OUT"
else
  OUT="$OUTPUT/public"
  python scripts/05_call_public.py --out "$OUT"
fi

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUT" "$REPO" submit_public.sbatch "$HSP65_AMPLICONS" "$RUNS/public_*/assembly/results"
