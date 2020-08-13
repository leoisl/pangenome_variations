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
from src.ConsistentPangenomeVariations import ConsistentPangenomeVariations


# setup
allele_mphf_filepath = snakemake.input.allele_mphf_filepath
consistent_pangenome_variations_filepath = snakemake.input.consistent_pangenome_variations
snps_dfs_filepaths = snakemake.input.snps_dfs_filepaths
deduplicated_snps_dfs_filepaths = snakemake.output.deduplicated_snps_dfs_filepaths
deduplicated_snps_dfs_text_filepaths = snakemake.output.deduplicated_snps_dfs_text_filepaths


# API usage
logging.info("Loading data...")
allele_mphf = AlleleMPHF.load(allele_mphf_filepath)
consistent_pangenome_variations = ConsistentPangenomeVariations.load(consistent_pangenome_variations_filepath)

logging.info("Converting to the deduplicated and enriched variations...")
for snps_df_filepath, deduplicated_snps_df_filepath, deduplicated_snps_df_text_filepath \
        in zip(snps_dfs_filepaths, deduplicated_snps_dfs_filepaths, deduplicated_snps_dfs_text_filepaths):
    deduplicated_snps_df = consistent_pangenome_variations.build_DeduplicatedVariationsDataframe_from_VarifierDataframe(
        snps_df_filepath, allele_mphf)
    deduplicated_snps_df.to_pickle(deduplicated_snps_df_filepath)
    deduplicated_snps_df.to_csv(deduplicated_snps_df_text_filepath, index=False)

logging.info("Done")