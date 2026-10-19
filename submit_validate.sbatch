#!/usr/bin/env bash
#SBATCH --job-name=hsp65-refalign
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/hsp65_refalign/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/hsp65-refalign"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa

cd "$REPO_ROOT"
python 01_validate_hsp65_refalign.py
