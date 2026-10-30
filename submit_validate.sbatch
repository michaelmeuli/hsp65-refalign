#!/usr/bin/env bash
# 01: Sanger hsp65 reads vs whole-genome species -> output/hsp65-refalign/results/ (reads.tsv, isolates.tsv, SUMMARY.txt)
#SBATCH --job-name=hsp65-refalign
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/hsp65-refalign/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/hsp65-refalign
OUTPUT=$OUT_HSP65/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

cd "$REPO"
python scripts/01_validate_hsp65_refalign.py --out "$OUTPUT"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_validate.sbatch "$SCREENING_MAP_RESULTS" "$SANGER_LINK_DIR" "$GTDB_REPRESENTATIVES"
