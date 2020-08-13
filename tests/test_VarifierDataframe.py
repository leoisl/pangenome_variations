from unittest import TestCase
from src.VarifierDataframe import VarifierDataframe
from io import StringIO
import pandas as pd


class TestVarifierDataframe(TestCase):
    ref_length = 20
    ref_pos = 3
    query_length = 10
    query_pos = 8
    rev_query_pos = 3

    @staticmethod
    def create_expected_Varifier_dataframe(
            query_strand_is_reversed: bool
    ) -> VarifierDataframe:
        if query_strand_is_reversed:
            query_pos = TestVarifierDataframe.rev_query_pos
        else:
            query_pos = TestVarifierDataframe.query_pos

        return VarifierDataframe(
            {
                "ref_chrom": ["ref_chrom"],
                "ref_pos": [TestVarifierDataframe.ref_pos],
                "ref_allele": ["G"],
                "ref_strand": [1],
                "ref_probe": ["AAAGTTT"],
                "ref_probe_interval": ["[3,4)"],
                "ref_len": [TestVarifierDataframe.ref_length],
                "query_chrom": ["query_chrom"],
                "query_pos": [query_pos],
                "query_allele": ["C"],
                "query_strand": [1],
                "query_probe": ["AAACTTT"],
                "query_probe_interval": ["[3,4)"],
                "query_len": [TestVarifierDataframe.query_length],
            }
        )

    @staticmethod
    def create_VCF_line_as_StringIO(
            query_strand_is_reversed: bool
    ) -> StringIO:
        query_strand = "-" if query_strand_is_reversed else "+"
        return StringIO(f"""#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
ref_chrom	{TestVarifierDataframe.ref_pos}	0	G	C	.	PASS	LENGTH_QRY={TestVarifierDataframe.query_length};LENGTH_REF={TestVarifierDataframe.ref_length};QNAME=query_chrom;QSTART={TestVarifierDataframe.query_pos};QSTRAND={query_strand}	GT:VFR_FILTER:VFR_REF_PROBE:VFR_REF_PROBE_ALLELE_INTERVAL:VFR_ALT_PROBE:VFR_ALT_PROBE_ALLELE_INTERVAL:VFR_ED_RA:VFR_ED_TR:VFR_ED_TA:VFR_ALLELE_LEN:VFR_ALLELE_MATCH_COUNT:VFR_ALLELE_MATCH_FRAC:VFR_IN_MASK:VFR_RESULT	1/1:PASS:AAAGTTT:[3,4):AAACTTT:[3,4):1:1:0:1:1:1.0:0:TP
""")

    def test_translateToFWDStrand_queryStrandIsForwardReturnsPositionsWithNoChange(
            self
    ):
        vcf_filehandler = TestVarifierDataframe.create_VCF_line_as_StringIO(query_strand_is_reversed=False)
        df = VarifierDataframe.load_from_varifier_VCF_core(vcf_filehandler)

        actual = df.translate_to_FWD_strand()
        expected = TestVarifierDataframe.create_expected_Varifier_dataframe(
            query_strand_is_reversed=False
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsReverseReturnsPositionsWithQueryChanged(
            self
    ):
        vcf_filehandler = TestVarifierDataframe.create_VCF_line_as_StringIO(query_strand_is_reversed=True)
        df = VarifierDataframe.load_from_varifier_VCF_core(vcf_filehandler)

        actual = df.translate_to_FWD_strand()
        expected = TestVarifierDataframe.create_expected_Varifier_dataframe(
            query_strand_is_reversed=True
        )

        assert actual.equals(expected)

    def test___get_ref_and_query_from_VarifierDataframe_filepath___absolute_path(self):
        filepath = "/home/leandro/git/pandora1_paper/analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = VarifierDataframe.get_ref_and_query_from_VarifierDataframe_filepath(filepath)
        assert ref == "CFT073" and query == "H131800734"

    def test___get_ref_and_query_from_VarifierDataframe_filepath___relative_path(self):
        filepath = "analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = VarifierDataframe.get_ref_and_query_from_VarifierDataframe_filepath(filepath)
        assert ref == "CFT073" and query == "H131800734"

    def test___get_ref_and_query_from_VarifierDataframe_filepath___local_path(self):
        filepath = "CFT073_and_H131800734.snps_df.pickle"
        ref, query = VarifierDataframe.get_ref_and_query_from_VarifierDataframe_filepath(filepath)
        assert ref == "CFT073" and query == "H131800734"

    def test___get_ref_and_query_from_VarifierDataframe_filepath___local_path___trivial_names(self):
        filepath = "A_and_B.snps_df.pickle"
        ref, query = VarifierDataframe.get_ref_and_query_from_VarifierDataframe_filepath(filepath)
        assert ref == "A" and query == "B"

    info = "LENGTH_QRY=4725092;LENGTH_REF=5155066;QNAME=1_pilon_pilon_pilon_pilon_pilon_pilon;QSTART=273255;QSTRAND=+;"
    def test___parse_field_from_info___LENGTH_QRY(self):
        expected = 4725092
        actual = VarifierDataframe.parse_field_from_info(TestVarifierDataframe.info, "LENGTH_QRY", int)
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___LENGTH_REF(self):
        expected = 5155066
        actual = VarifierDataframe.parse_field_from_info(TestVarifierDataframe.info, "LENGTH_REF", int)
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___QNAME(self):
        expected = "1_pilon_pilon_pilon_pilon_pilon_pilon"
        actual = VarifierDataframe.parse_field_from_info(TestVarifierDataframe.info, "QNAME")
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___QSTART(self):
        expected = 273255
        actual = VarifierDataframe.parse_field_from_info(TestVarifierDataframe.info, "QSTART", int)
        self.assertEqual(expected, actual)

    def test___parse_field_from_info___QSTRAND(self):
        expected = "+"
        actual = VarifierDataframe.parse_field_from_info(TestVarifierDataframe.info, "QSTRAND")
        self.assertEqual(expected, actual)

    def test___load_from_varifier_VCF_core(self):
        vcf_string_io = StringIO(
"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=0,length=4824417>
##contig=<ID=1,length=173787>
##contig=<ID=2,length=62174>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##bcftools_normVersion=1.10.2 (pysam)+htslib-1.10.2 (pysam)
##bcftools_normCommand=norm -c x -d any -f Escherichia_coli_MSB1_1A.ref.fa -o /hps/nobackup2/iqbal/leandro/temp/tmpa_qmdli3 test_out/03.probe_filtered.vcf; Date=Wed Aug 12 15:07:59 2020
##bcftools_normCommand=norm -c x -d any -f Escherichia_coli_MSB1_1A.ref.fa -o /hps/nobackup2/iqbal/leandro/temp/tmpxver718t test_out/04.truth.vcf; Date=Wed Aug 12 15:08:01 2020
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
ref_chrom1       345     0       A       G       .       PASS    LENGTH_QRY=4905296;LENGTH_REF=4824417;QNAME=query_chrom_0;QSTART=1060;QSTRAND=-     GT:VFR_FILTER:VFR_REF_PROBE:VFR_REF_PROBE_ALLELE_INTERVAL:VFR_ALT_PROBE:VFR_ALT_PROBE_ALLELE_INTERVAL:VFR_ED_RA:VFR_ED_TR:VFR_ED_TA:VFR_ALLELE_LEN:VFR_ALLELE_MATCH_COUNT:VFR_ALLELE_MATCH_FRAC:VFR_IN_MASK:VFR_RESULT  1/1:PASS:CGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAAC:[100,101):CGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCGGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAAC:[100,101):1:1:0:1:1:1.0:0:TP
ref_chrom2       366     1       C       T       .       PASS    LENGTH_QRY=4905296;LENGTH_REF=4824417;QNAME=query_chrom_1;QSTART=1039;QSTRAND=+     GT:VFR_FILTER:VFR_REF_PROBE:VFR_REF_PROBE_ALLELE_INTERVAL:VFR_ALT_PROBE:VFR_ALT_PROBE_ALLELE_INTERVAL:VFR_ED_RA:VFR_ED_TR:VFR_ED_TA:VFR_ALLELE_LEN:VFR_ALLELE_MATCH_COUNT:VFR_ALLELE_MATCH_FRAC:VFR_IN_MASK:VFR_RESULT  1/1:PASS:CGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCC:[100,101):CGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCTCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCC:[100,101):1:1:0:1:1:1.0:0:TP
""")
        expected_csv_string_io = StringIO(
"""ref_chrom,ref_pos,ref_allele,ref_strand,ref_probe,ref_probe_interval,ref_len,query_chrom,query_pos,query_allele,query_strand,query_probe,query_probe_interval,query_len
ref_chrom1,345,A,1,CGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAAC,"[100,101)",4824417,query_chrom_0,4904237,G,1,CGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCGGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAAC,"[100,101)",4905296
ref_chrom2,366,C,1,CGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCC,"[100,101)",4824417,query_chrom_1,1039,T,1,CGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCTCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCC,"[100,101)",4905296
""")
        actual_df = VarifierDataframe.load_from_varifier_VCF_core(vcf_string_io)
        expected_df = VarifierDataframe(pd.read_csv(expected_csv_string_io))
        self.assertTrue(actual_df.equals(expected_df))
