import itertools
import pandas as pd
from snakemake.utils import min_version

min_version("5.1.0")

# ======================================================
# Config files and vars
# ======================================================
configfile: "config.yaml"
output_folder = config['output_folder']
samples_file = config["samples"]
flank_length = config["flank_length"]

samples = pd.read_csv(samples_file)
samples.rename(columns={"reference": "reference_assembly"}, inplace=True)
samples = samples.set_index(["sample_id"], drop=False)
sample_pairs = [(sample1, sample2) for sample1, sample2 in itertools.combinations(sorted(samples["sample_id"]), r=2)]
sample_pairs_as_str = [f"{sample1}/{sample1}_and_{sample2}" for sample1, sample2 in sample_pairs]



# ======================================================
# Rules
# ======================================================
rule all:
    input:
        expand(output_folder + "/truth_probesets/{sample_pairs_as_str}.truth_probeset.fa", sample_pairs_as_str=sample_pairs_as_str)

rule make_pairwise_snps_df:
    input:
         truth1=lambda wildcards: samples.xs(wildcards.sample1)["reference_assembly"],
         truth2=lambda wildcards: samples.xs(wildcards.sample2)["reference_assembly"],
    output:
          snps_df=output_folder + "/snps_dfs/{sample1}/{sample1}_and_{sample2}.snps_df.pickle",
          snps_df_text=output_folder + "/snps_dfs/{sample1}/{sample1}_and_{sample2}.snps_df.csv",
          aligned_bases_percentage_sample_1=output_folder + "/dnadiff_reports/{sample1}/{sample1}_and_{sample2}.aligned_bases_percentage",
          aligned_bases_percentage_sample_2=output_folder + "/dnadiff_reports/{sample2}/{sample1}_and_{sample2}.aligned_bases_percentage"
    params:
          flank_length=flank_length
    shadow:
          "shallow"
    threads: 1
    resources:
             mem_mb=lambda wildcards, attempt: 4000 * attempt
    log:
       "logs/make_pairwise_snps_df/{sample1}_and_{sample2}.log"
    script:
          "scripts/make_pairwise_snps_df.py"


rule make_allele_mphf:
    input:
         snps_dfs=expand(output_folder + "/snps_dfs/{sample_pair_as_str}.snps_df.pickle",
                         sample_pair_as_str=sample_pairs_as_str)
    output:
          number_of_alleles=output_folder + "/allele_mphf/number_of_alleles.txt",
          alelle_mphf=output_folder + "/allele_mphf/allele_mphf.pickle",
    threads: 1
    resources:
             mem_mb=lambda wildcards, attempt: 16000 * attempt
    log:
       "logs/make_allele_mphf.log"
    script:
          "scripts/make_allele_mphf.py"


rule make_pairwise_variants_mphf:
    input:
         snps_dfs=expand(output_folder + "/snps_dfs/{sample_pair_as_str}.snps_df.pickle",
                         sample_pair_as_str=sample_pairs_as_str),
         allele_mphf=rules.make_allele_mphf.output.alelle_mphf
    output:
          number_of_pairwise_variants=output_folder + "/pairwise_variants_mphf/number_of_pairwise_variants.txt",
          pairwise_variants_mphf=output_folder + "/pairwise_variants_mphf/pairwise_variants_mphf.pickle",
          pairwise_variation_id_to_alleles_id=output_folder + "/pairwise_variants_mphf/pairwise_variation_id_to_alleles_id.pickle",
    threads: 1
    resources:
             mem_mb=lambda wildcards, attempt: 16000 * attempt
    log:
       "logs/make_pairwise_variants_mphf.log"
    script:
          "scripts/make_pairwise_variants_mphf.py"


rule deduplicate_pairwise_snps:
    input:
         number_of_alleles_filepath=rules.make_allele_mphf.output.number_of_alleles,
         pairwise_variation_id_to_alleles_id_filepath=rules.make_pairwise_variants_mphf.output.pairwise_variation_id_to_alleles_id,
    output:
          pangenome_variations_defined_by_allele_ids=output_folder + "/deduplicate_pairwise_snps/pangenome_variations_defined_by_allele_ids.pickle",
    threads: 1
    resources:
             mem_mb=lambda wildcards, attempt: 16000 * attempt
    log:
       "logs/deduplicate_pairwise_snps.log"
    script:
          "scripts/deduplicate_pairwise_snps.py"


rule convert_pangenome_variations_to_deduplicated_snps_df:
    input:
         pangenome_variations_defined_by_allele_ids_filepath=rules.deduplicate_pairwise_snps.output.pangenome_variations_defined_by_allele_ids,
         allele_mphf_filepath=rules.make_allele_mphf.output.alelle_mphf,
         snps_dfs_filepaths=expand(output_folder + "/snps_dfs/{sample_pair_as_str}.snps_df.pickle",
                                   sample_pair_as_str=sample_pairs_as_str),
    output:
          deduplicated_snps_dfs_filepaths=expand(
              output_folder + "/deduplicated_snps_dfs/{sample_pair_as_str}.snps_df.pickle",
              sample_pair_as_str=sample_pairs_as_str),
          deduplicated_snps_dfs_text_filepaths=expand(
              output_folder + "/deduplicated_snps_dfs/{sample_pair_as_str}.snps_df.csv",
              sample_pair_as_str=sample_pairs_as_str)
    threads: 1
    resources:
             mem_mb=lambda wildcards, attempt: 16000 * attempt
    log:
       "logs/convert_pangenome_variations_to_deduplicated_snps_df.log"
    script:
          "scripts/convert_pangenome_variations_to_deduplicated_snps_df.py"


rule make_truth_probeset:
    input:
         deduplicated_snps_df=output_folder + "/deduplicated_snps_dfs/{sample1}/{sample1}_and_{sample2}.snps_df.pickle",
    output:
          probeset1=output_folder + "/truth_probesets/{sample1}/{sample1}_and_{sample2}.truth_probeset.fa",
          probeset2=output_folder + "/truth_probesets/{sample2}/{sample1}_and_{sample2}.truth_probeset.fa",
    threads: 1
    resources:
             mem_mb=lambda wildcards, attempt: 4000 * attempt
    log:
       "logs/make_truth_probeset/{sample1}_and_{sample2}.log"
    script:
          "scripts/make_truth_probeset.py"
