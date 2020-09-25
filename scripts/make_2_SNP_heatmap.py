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
two_SNP_heatmap_csv_filename = snakemake.output.two_SNP_heatmap_csv


# API usage
logging.info("Loading data...")
consistent_pangenome_variations = ConsistentPangenomeVariations.load(consistent_pangenome_variations_filename)

logging.info("Building csv...")
consistent_pangenome_variations_in_two_samples_only = [variation
                                                       for variation in
                                                       consistent_pangenome_variations.consistent_pangenome_variations
                                                       if variation.get_number_of_samples()==2]

ids = []
first_samples = []
second_samples = []

for variation in consistent_pangenome_variations_in_two_samples_only:
    samples = sorted(list(variation.get_set_of_samples()))
    ids.append(variation.id)
    first_samples.append(samples[0])
    second_samples.append(samples[1])

two_SNP_heatmap_df = pd.DataFrame(data={
    "PANGENOME_VARIATION_ID": ids,
    "FIRST_SAMPLE": first_samples,
    "SECOND_SAMPLE": second_samples
})


# output
logging.info("Saving csv...")
two_SNP_heatmap_df.to_csv(two_SNP_heatmap_csv_filename, index=False)

logging.info("Done")