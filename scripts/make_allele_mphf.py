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
snps_dfs = snakemake.input.snps_dfs
number_of_alleles = snakemake.output.number_of_alleles
alelle_mphf = snakemake.output.alelle_mphf

allele_mphf = AlleleMPHF.build_from_list_of_snps_dfs_filepaths(snps_dfs)
with open(number_of_alleles, "w") as number_of_alleles_fh, \
        open(alelle_mphf, "wb") as allele_mphf_fh:
    allele_mphf.dump(number_of_alleles_fh, allele_mphf_fh)
