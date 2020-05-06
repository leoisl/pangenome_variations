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
from io import StringIO
from src.mummer import ShowSnps, Nucmer, DeltaFilter, NucmerError, GetReportFromDeltaFile


def generate_mummer_snps(
        reference: Path,
        query: Path,
        prefix: Path = Path("out"),
        flank_width: int = 0,
        indels: bool = True,
        print_header: bool = True,
) -> StringIO:
    nucmer_params = "--maxmatch"
    nucmer = Nucmer(reference, query, str(prefix), extra_params=nucmer_params)
    nucmer_result = nucmer.run()
    if nucmer_result.returncode != 0:
        raise NucmerError(nucmer_result.stderr.decode())

    deltafile = Path(str(prefix) + ".delta")
    deltafilter_params = "-1"
    deltafilter = DeltaFilter(deltafile, extra_params=deltafilter_params)
    deltafilter_result = deltafilter.run()
    deltafilter_result.check_returncode()

    filtered_deltafile = prefix.with_suffix(".delta1")
    _ = filtered_deltafile.write_text(deltafilter_result.stdout.decode())

    showsnps_params = "-rlTC"
    showsnps = ShowSnps(
        filtered_deltafile,
        context=flank_width,
        extra_params=showsnps_params,
        indels=indels,
        print_header=print_header,
    )
    showsnps_result = showsnps.run()
    showsnps_result.check_returncode()
    showsnps_content = showsnps_result.stdout.decode()

    snpsfile = prefix.with_suffix(".snps")
    _ = snpsfile.write_text(showsnps_content)

    return StringIO(showsnps_content)


def get_dnadiff_report(
        prefix: Path = Path("out")
) -> str:
    deltafile = Path(str(prefix) + ".delta")
    get_report_from_delta_file = GetReportFromDeltaFile(deltafile)
    report_as_str = get_report_from_delta_file.get_report()
    return report_as_str


# ==================================================================================
# MAIN
# ==================================================================================

# setup
query1: Path = Path(snakemake.input.truth1)
query2: Path = Path(snakemake.input.truth2)
query1_name: str = snakemake.wildcards.sample1
query2_name: str = snakemake.wildcards.sample2
prefix: Path = Path(f"{query1_name}_and_{query2_name}")
flank_width: int = snakemake.params.flank_length


# API usage
logging.info("Generating mummer snps")
mummer_snps: StringIO = generate_mummer_snps(
    reference=query1, query=query2, prefix=prefix, flank_width=flank_width, indels=False
)

logging.info("Generating dnadiff report")
dna_diff_report_as_str = get_dnadiff_report(prefix)
ref_aligned_bases_percentage, query_aligned_bases_percentage = \
    GetReportFromDeltaFile.get_ref_and_query_aligned_bases_percentage(dna_diff_report_as_str)
aligned_bases_percentage_sample_1 = Path(snakemake.output.aligned_bases_percentage_sample_1)
aligned_bases_percentage_sample_2 = Path(snakemake.output.aligned_bases_percentage_sample_2)
aligned_bases_percentage_sample_1.write_text(str(ref_aligned_bases_percentage))
aligned_bases_percentage_sample_2.write_text(str(query_aligned_bases_percentage))

logging.info("Converting show-snps output to dataframe")
snps_df = ShowSnps.to_dataframe(mummer_snps)
snps_df = snps_df.translate_to_FWD_strand()
snps_df.to_pickle(snakemake.output.snps_df)
snps_df.to_csv(snakemake.output.snps_df_text, index=False)
logging.info(f"Done")
