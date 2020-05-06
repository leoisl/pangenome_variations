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
snps_dfs_filepaths = snakemake.input.snps_dfs_filepaths
deduplicated_snps_dfs_filepaths = snakemake.output.deduplicated_snps_dfs_filepaths
deduplicated_snps_dfs_text_filepaths = snakemake.output.deduplicated_snps_dfs_text_filepaths


# API usage
allele_mphf = AlleleMPHF.load(allele_mphf_filepath)
with open(pangenome_variations_defined_by_allele_ids_filepath, "rb") as pangenome_variations_defined_by_allele_ids_file:
    pangenome_variations_defined_by_allele_ids = pickle.load(pangenome_variations_defined_by_allele_ids_file)

pangenome_variations = PangenomeVariations.build_from_pangenome_variations_defined_by_allele_ids(
    pangenome_variations_defined_by_allele_ids, allele_mphf
)
consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations, filter_for_biallelic=True)
logging.info(f"Number of pangenome variations: {consistent_pangenome_variations.number_of_pangenome_variations}")
logging.info(f"Number of consistent pangenome variations: {consistent_pangenome_variations.number_of_consistent_pangenome_variations}")
logging.info(f"Number of consistent biallelic pangenome variations: {consistent_pangenome_variations.number_of_consistent_biallelic_pangenome_variations}")

# write the enriched and filtered deduplicated variations
for snps_df_filepath, deduplicated_snps_df_filepath, deduplicated_snps_df_text_filepath \
        in zip(snps_dfs_filepaths, deduplicated_snps_dfs_filepaths, deduplicated_snps_dfs_text_filepaths):
    deduplicated_snps_df = consistent_pangenome_variations.build_DeduplicatedVariationsDataframe_from_ShowSNPsDataframe(
        snps_df_filepath, allele_mphf)
    deduplicated_snps_df.to_pickle(deduplicated_snps_df_filepath)
    deduplicated_snps_df.to_csv(deduplicated_snps_df_text_filepath, index=False)
