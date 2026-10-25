#!/usr/bin/env bash
# sbatch submit_public.sbatch        GTDB-kansasii samples of the public BioProjects -> output/hsp65_refalign/public/
# sbatch submit_public.sbatch ctl    positive control: first 3 GTDB-persicum samples per project -> output/hsp65_refalign/public_ctl/
#SBATCH --job-name=hsp65-public
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/hsp65_refalign/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/hsp65-refalign"
OUT_ROOT="/shares/sander.imm.uzh/MM/kansasii/output/hsp65_refalign"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

cd "$REPO_ROOT"
if [ "${1:-}" = "ctl" ]; then
  OUT="$OUT_ROOT/public_ctl"
  python 05_call_public.py --limit 3 --species "s__Mycobacterium persicum" --out "$OUT"
else
  OUT="$OUT_ROOT/public"
  python 05_call_public.py --out "$OUT"
fi

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "$OUT" \
  /shares/sander.imm.uzh/MM/kansasii/repos/hsp65-refalign submit_public.sbatch \
  /shares/sander.imm.uzh/MM/kansasii/data/gtdb_genomes/Mycobacteriaceae/kansasii_complex_hsp65_amplicons
