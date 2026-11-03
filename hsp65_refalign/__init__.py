"""Repo-specific constants; project-wide paths, species and colours are in `kansasii_paths` (shared, vendored)."""
from .kansasii_paths import OUT_HSP65

PUBLIC_MAX_CONTAMINATION = 5.0  # CheckM contamination (%) above which a public sample is left out of every result
OUTPUT = OUT_HSP65 / "results"
OUTPUT_PUBLIC = OUTPUT / "public"  # 05: public BioProject assemblies
