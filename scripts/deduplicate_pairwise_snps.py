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
from src.DeduplicationGraph import DeduplicationGraph
import pickle


# setup
number_of_alleles_filepath = Path(snakemake.input.number_of_alleles_filepath)
pairwise_variation_id_to_alleles_id_filepath = snakemake.input.pairwise_variation_id_to_alleles_id_filepath
pangenome_variations_defined_by_allele_ids_filepath = snakemake.output.pangenome_variations_defined_by_allele_ids

number_of_alleles = int(number_of_alleles_filepath.read_text().strip())
logging.info("Loading pairwise_variation_id_to_alleles_id...")
with open(pairwise_variation_id_to_alleles_id_filepath, "rb") as pairwise_variation_id_to_alleles_id_file:
    pairwise_variation_id_to_alleles_id = pickle.load(pairwise_variation_id_to_alleles_id_file)

# API usage
logging.info("Creating DeduplicationGraph...")
deduplication_graph = DeduplicationGraph(number_of_alleles, pairwise_variation_id_to_alleles_id)
logging.info("Getting pangenome variations defined by allele ids...")
pangenome_variations = deduplication_graph.get_pangenome_variations_defined_by_allele_ids()
logging.info("Saving pangenome variations defined by allele ids...")
with open(pangenome_variations_defined_by_allele_ids_filepath, "wb") as pangenome_variations_defined_by_allele_ids_file:
    pickle.dump(pangenome_variations, pangenome_variations_defined_by_allele_ids_file)
logging.info("Done")