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
from src.AlleleMPHF import AlleleMPHF


# setup
snps_dfs_filepaths = snakemake.input.snps_dfs
number_of_alleles_filepath = snakemake.output.number_of_alleles
alelle_mphf_filepath = snakemake.output.alelle_mphf


# API usage
allele_mphf = AlleleMPHF.build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths)
logging.info("Saving AlleleMPHF...")
allele_mphf.save(number_of_alleles_filepath, alelle_mphf_filepath)
logging.info("Done!")
