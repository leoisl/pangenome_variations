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
import pandas as pd


# setup
deduplicated_snps_df_filepath = snakemake.input.deduplicated_snps_df
probeset1_filepath = snakemake.output.probeset1
probeset2_filepath = snakemake.output.probeset2


# API usage
logging.info("Reading show-snps dataframe")
snps_df = pd.read_pickle(deduplicated_snps_df_filepath)

logging.info("Creating probes from dataframe")
query1_truth_probes, query2_truth_probes = snps_df.get_probes()

logging.info("Writing output files")
query1_truth_probes_path: Path = Path(probeset1_filepath)
query2_truth_probes_path: Path = Path(probeset2_filepath)
query1_truth_probes_path.write_text(query1_truth_probes)
query2_truth_probes_path.write_text(query2_truth_probes)

logging.info(f"Done")
