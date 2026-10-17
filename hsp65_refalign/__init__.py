import os
from pathlib import Path

SPECIES = ["kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"]

# Mkan329-126 to -133 are reference strains, not patient isolates; they are excluded from every analysis
REFERENCE_STRAINS = frozenset(f"Mkan329-{i:03d}" for i in range(126, 134))

# fixed colour per species (same as mlsa-kansasii) so figures stay comparable
SPECIES_COLORS = {"kansasii": "#1b9e77", "persicum": "#d95f02", "pseudokansasii": "#7570b3", "innocens": "#e7298a",
                  "attenuatum": "#66a61e", "ostraviense": "#e6ab02", "gastri": "#a6761d"}

# Everything lives under the shared kansasii tree; repos/<this repo> sits two levels below it.
ROOT = Path(os.environ.get("KANSASII_ROOT", Path(__file__).resolve().parent.parent.parent.parent))
DATA = ROOT / "data"
SANGER_LINK_DIR = DATA / "sanger" / "seq_kansasii"
# one GTDB representative genome per species (the references of the plain read alignment)
GTDB_REPRESENTATIVES = DATA / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_gtdb_representatives"
# in-silico hsp65 amplicons (441 bp) of the 72 public complex genomes, <accession>_hsp65.fasta
HSP65_AMPLICONS = DATA / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_hsp65_amplicons"
SCREENING_RESULTS = ROOT / "output" / "screening_map_results.csv"
ASSEMBLIES = ROOT / "runs" / "mkan329" / "assembly" / "results"
OUTPUT = ROOT / "output" / "hsp65_refalign"
