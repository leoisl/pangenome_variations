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
from src.ConsistentPangenomeVariations import ConsistentPangenomeVariations
import pandas as pd


# setup
consistent_pangenome_variations_filename = snakemake.input.consistent_pangenome_variations
id_and_number_of_samples_csv_filename = snakemake.output.id_and_number_of_samples_csv


# API usage
logging.info("Loading data...")
consistent_pangenome_variations = ConsistentPangenomeVariations.load(consistent_pangenome_variations_filename)

logging.info("Building csv...")
ids = [variation.id for variation in consistent_pangenome_variations.consistent_pangenome_variations]
numbers_of_samples = [variation.get_number_of_samples()
                      for variation in consistent_pangenome_variations.consistent_pangenome_variations]
id_and_number_of_samples_df = pd.DataFrame(data={
    "PVID": ids,
    "NUMBER_OF_SAMPLES": numbers_of_samples})


# output
logging.info("Saving csv...")
id_and_number_of_samples_df.to_csv(id_and_number_of_samples_csv_filename, index=False)

logging.info("Done")