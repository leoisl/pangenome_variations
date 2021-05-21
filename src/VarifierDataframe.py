import pandas as pd
from typing import Tuple, TextIO
from pathlib import Path
import re


class VarifierDataframe(pd.DataFrame):
    @property
    def _constructor(self):
        return VarifierDataframe

    def translate_to_FWD_strand(self) -> "VarifierDataframe":
        def fix_position(position: int, strand_aln: int, length: int) -> int:
            if strand_aln == 1:
                return position
            else:
                return length - position + 1

        def translate_to_FWD_strand_core(line: pd.Series) -> pd.Series:
            # ref is always in the FWD strand, there is no need to translate it
            line.query_pos = fix_position(
                line.query_pos, line.query_strand, line.query_len
            )
            line.query_strand = 1
            return line

        return self.apply(translate_to_FWD_strand_core, axis=1)

    @staticmethod
    def parse_field_from_info(info: str, field: str, type=str):
        field_matcher = re.compile(f".*{field}=(.+?);")
        result = field_matcher.match(info).group(1)
        return type(result)

    @staticmethod
    def get_ref_and_query_from_VarifierDataframe_filepath(VarifierDataframe_filepath: str) -> Tuple[str, str]:
        VarifierDataframe_filepath = Path(VarifierDataframe_filepath)
        VarifierDataframe_filename = VarifierDataframe_filepath.name
        matches = re.match(r"(.*)_and_(.*).snps_df.pickle", VarifierDataframe_filename)
        ref = matches.group(1)
        query = matches.group(2)
        return ref, query


    # Note: not tested
    @staticmethod
    def load_pickled(df_filepath: str) -> "VarifierDataframe":
        varifier_df = pd.read_pickle(df_filepath)
        varifier_df = varifier_df.translate_to_FWD_strand()
        return varifier_df

    # Note: not tested
    @staticmethod
    def load_from_varifier_VCF(vcf_filepath: str) -> "VarifierDataframe":
        with open(vcf_filepath) as vcf_filehandler:
            return VarifierDataframe.load_from_varifier_VCF_core(vcf_filehandler)

    @staticmethod
    def load_from_varifier_VCF_core(vcf_filehandler: TextIO) -> "VarifierDataframe":
        """
        Example of VCF line:
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
        0	28	0	A	T	.	PASS	LENGTH_QRY=4824417;LENGTH_REF=4905296;QNAME=0;QSTART=1377;QSTRAND=-	GT:VFR_FILTER:VFR_REF_PROBE:VFR_ALT_PROBE:VFR_ED_RA:VFR_ED_TR:VFR_ED_TA:VFR_ALLELE_LEN:VFR_ALLELE_MATCH_COUNT:VFR_ALLELE_MATCH_FRAC:VFR_IN_MASK:VFR_RESULT	1/1:PASS:TTACGATGACAATGTTCTGATTAAATTAGAAAAATCTTCTTTGATATCGTGGCTCTCTTCACGCAACTGCTCGATCTTCCGGCAGGCATGAAGCACCGTCGTGTGGTCACGACCACCAAAGGCATCGC:TTACGATGACAATGTTCTGATTAAATTTGAAAAATCTTCTTTGATATCGTGGCTCTCTTCACGCAACTGCTCGATCTTCCGGCAGGCATGAAGCACCGTCGTGTGGTCACGACCACCAAAGGCATCGC:1:1:0:1:1:1.0:0:TP

        """
        ref_chrom = []
        ref_pos = []
        ref_allele = []
        ref_strand = []
        ref_probe = []
        ref_probe_interval = []
        ref_len = []

        query_chrom = []
        query_pos = []
        query_allele = []
        query_strand = []
        query_probe = []
        query_probe_interval = []
        query_len = []

        for line in vcf_filehandler:
            line = line.strip()
            is_header = line.startswith("#")
            if not is_header:
                line_split = line.split()

                ref_chrom.append(line_split[0])
                ref_pos.append(int(line_split[1]))
                ref_allele.append(line_split[3])
                ref_strand.append(1)
                query_allele.append(line_split[4])

                info = line_split[7]
                info += ";" # makes regex matching a bit easier
                query_len.append(VarifierDataframe.parse_field_from_info(info, "LENGTH_QRY", int))
                ref_len.append(VarifierDataframe.parse_field_from_info(info, "LENGTH_REF", int))
                query_chrom.append(VarifierDataframe.parse_field_from_info(info, "QNAME"))
                query_pos.append(VarifierDataframe.parse_field_from_info(info, "QSTART", int))
                binary_query_strand = 1 if VarifierDataframe.parse_field_from_info(info, "QSTRAND") == "+" else -1
                query_strand.append(binary_query_strand)

                sample_data = line_split[9].split(":")
                ref_probe.append(sample_data[2])
                ref_probe_interval.append(sample_data[3])
                query_probe.append(sample_data[4])
                query_probe_interval.append(sample_data[5])

        varifier_df = VarifierDataframe(data={
            "ref_chrom": ref_chrom,
            "ref_pos": ref_pos,
            "ref_allele": ref_allele,
            "ref_strand": ref_strand,
            "ref_probe": ref_probe,
            "ref_probe_interval": ref_probe_interval,
            "ref_len": ref_len,
            "ref_original_pos": ref_pos,
            "ref_original_strand": ref_strand,
            "query_chrom": query_chrom,
            "query_pos": query_pos,
            "query_allele": query_allele,
            "query_strand": query_strand,
            "query_probe": query_probe,
            "query_probe_interval": query_probe_interval,
            "query_len": query_len,
            "query_original_pos": query_pos,
            "query_original_strand": query_strand,
        })

        varifier_df = varifier_df.translate_to_FWD_strand()
        return varifier_df

