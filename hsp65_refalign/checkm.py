"""CheckM contamination of the public samples (same rule as gyrA-refalign: above PUBLIC_MAX_CONTAMINATION a sample is left out)."""
from __future__ import annotations

import re
from pathlib import Path


def checkm_contamination(assembly_dir: Path, sample: str) -> float | None:
    """CheckM contamination (%) of a public sample from its 3_quality/checkM table; None if the table is missing."""
    f = assembly_dir / sample / "3_quality" / "checkM" / "checkm_output.tsv"
    if not f.exists():
        return None
    for line in f.read_text().splitlines():
        if re.match(rf"\s+{re.escape(sample)}\s", line):
            return float(line.split()[-2])  # ... Completeness  Contamination  Strain heterogeneity
    return None
