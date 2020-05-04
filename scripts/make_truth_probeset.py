import sys
from pathlib import Path

sys.path.append(str(Path().absolute()))
import logging

log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)
import pickle

# ==================================================================================
# MAIN
# ==================================================================================
logging.info("Reading show-snps dataframe")
with open(snakemake.input.deduplicated_snps_df, "rb") as snps_df_fh:
    snps_df = pickle.load(snps_df_fh)

logging.info("Creating probes from dataframe")
query1_truth_probes, query2_truth_probes = snps_df.get_probes()

logging.info("Writing output files")
query1_truth_probes_path: Path = Path(snakemake.output.probeset1)
query2_truth_probes_path: Path = Path(snakemake.output.probeset2)
query1_truth_probes_path.write_text(query1_truth_probes)
query2_truth_probes_path.write_text(query2_truth_probes)

logging.info(f"Done")
