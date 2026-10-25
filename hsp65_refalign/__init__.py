"""Repo-specific constants; project-wide paths, species and colours are in `kansasii_paths` (shared, vendored)."""
from .kansasii_paths import OUTPUT_ROOT

OUTPUT = OUTPUT_ROOT / "hsp65_refalign"
OUTPUT_PUBLIC = OUTPUT / "public"  # 05: public BioProject assemblies
