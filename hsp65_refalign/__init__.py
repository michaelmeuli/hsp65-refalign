import os
from pathlib import Path

SPECIES = ["kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"]

# Everything lives under the shared kansasii tree; repos/<this repo> sits two levels below it.
ROOT = Path(os.environ.get("KANSASII_ROOT", Path(__file__).resolve().parent.parent.parent.parent))
DATA = ROOT / "data"
SANGER_LINK_DIR = DATA / "sanger" / "seq_kansasii"
# one GTDB representative genome per species (the references of the plain read alignment)
GTDB_REPRESENTATIVES = DATA / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_gtdb_representatives"
SCREENING_RESULTS = ROOT / "output" / "screening_map_results.csv"
ASSEMBLIES = ROOT / "runs" / "mkan329" / "assembly" / "results"
OUTPUT = ROOT / "output" / "hsp65_refalign"
