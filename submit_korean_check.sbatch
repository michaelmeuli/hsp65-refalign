#!/usr/bin/env bash
# 04: Korean-type (persicum-like) hsp65 among the project isolates' own assemblies -> output/hsp65_refalign/korean_type_check.tsv
#SBATCH --job-name=hsp65-korean-check
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
python scripts/04_korean_type_check.py --out "$OUTPUT"

# provenance; MANIFEST_NAME keeps it apart from submit_validate's MANIFEST.txt in the same directory
MANIFEST_NAME=MANIFEST_korean_type_check.txt bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_korean_check.sbatch \
  "$SCREENING_MAP_RESULTS" "$ASSEMBLIES" "$HSP65_AMPLICONS"
