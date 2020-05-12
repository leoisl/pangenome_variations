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
from src.AlleleMPHF import AlleleMPHF
from src.PangenomeVariations import PangenomeVariations
from src.ConsistentPangenomeVariations import ConsistentPangenomeVariations


# setup
allele_mphf_filepath = snakemake.input.allele_mphf_filepath
pangenome_variations_defined_by_allele_ids_filepath = snakemake.input.pangenome_variations_defined_by_allele_ids_filepath
consistent_pangenome_variations_filepath = snakemake.output.consistent_pangenome_variations

# API usage
logging.info("Loading data...")
allele_mphf = AlleleMPHF.load(allele_mphf_filepath)
with open(pangenome_variations_defined_by_allele_ids_filepath, "rb") as pangenome_variations_defined_by_allele_ids_file:
    pangenome_variations_defined_by_allele_ids = pickle.load(pangenome_variations_defined_by_allele_ids_file)

logging.info("Building Pangenome Variations...")
pangenome_variations = PangenomeVariations.build_from_pangenome_variations_defined_by_allele_ids(
    pangenome_variations_defined_by_allele_ids, allele_mphf
)

logging.info("Building Consistent Biallelic Pangenome Variations...")
consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations, filter_for_biallelic=True)
logging.info(f"Number of pangenome variations: {consistent_pangenome_variations.number_of_pangenome_variations}")
logging.info(f"Number of consistent pangenome variations: {consistent_pangenome_variations.number_of_consistent_pangenome_variations}")
logging.info(f"Number of consistent biallelic pangenome variations: {consistent_pangenome_variations.number_of_consistent_biallelic_pangenome_variations}")

logging.info("Saving the Consistent Biallelic Pangenome Variations...")
consistent_pangenome_variations.save(consistent_pangenome_variations_filepath)

logging.info("Done")