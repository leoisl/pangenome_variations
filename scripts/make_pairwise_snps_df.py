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
from snakemake import shell
from src.VarifierDataframe import VarifierDataframe

# ==================================================================================
# MAIN
# ==================================================================================

# setup
query1: Path = Path(snakemake.input.truth1)
query2: Path = Path(snakemake.input.truth2)
query1_name: str = snakemake.wildcards.sample1
query2_name: str = snakemake.wildcards.sample2
prefix: Path = Path(f"{query1_name}_and_{query2_name}")
flank_length: int = snakemake.params.flank_length
varifier_vcf = snakemake.output.varifier_vcf
snps_df_filepath = snakemake.output.snps_df
snps_df_text_filepath = snakemake.output.snps_df_text


# API usage
logging.info("Generating varifier snps")
shell(f"varifier make_truth_vcf --snps_only --output_probes_in_VCF --detailed_VCF --flank_length {flank_length} {query2} {query1} {prefix}_varifier_output")
shell(f"cp {prefix}_varifier_output/04.truth.vcf {varifier_vcf}")

logging.info("Converting varifier VCF to dataframe")
snps_df = VarifierDataframe.load_from_varifier_VCF(varifier_vcf)
snps_df.to_pickle(snps_df_filepath)
snps_df.to_csv(snps_df_text_filepath, index=False)

logging.info(f"Done")
