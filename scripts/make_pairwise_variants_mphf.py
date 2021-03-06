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
from src.PairwiseVariationMPHF import PairwiseVariationMPHF


# setup
snps_dfs_filenames = snakemake.input.snps_dfs
allele_mphf_filename = snakemake.input.allele_mphf
number_of_pairwise_variants_filename = snakemake.output.number_of_pairwise_variants
pairwise_variants_mphf_filename = snakemake.output.pairwise_variants_mphf
pairwise_variation_id_to_alleles_id_filename = snakemake.output.pairwise_variation_id_to_alleles_id


# API usage
pairwise_variants_mphf = PairwiseVariationMPHF.build_from_list_of_snps_dfs_filepaths(snps_dfs_filenames,
                                                                                     allele_mphf_filename)
logging.info("Saving PairwiseVariationMPHF...")
pairwise_variants_mphf.save(number_of_pairwise_variants_filename, pairwise_variants_mphf_filename,
                            pairwise_variation_id_to_alleles_id_filename)
logging.info("Done!")
