import uuid
from io import StringIO

import pytest

from src.mummer import *

REF = Path("test_cases/ref.fa")
QUERY = Path("test_cases/query.fa")
DELTA = Path("test_cases/out.delta")
DELTA1 = Path("test_cases/out.delta1")


class TestNucmer:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        reference = Path("foo")
        query = Path("bar")
        with pytest.raises(NucmerError):
            nucmer = Nucmer(reference, query)

    def test_generateCommand_defaultArgsHasNoExtraParamsAndHasOutAsPrefix(self):
        reference = REF
        query = QUERY
        nucmer = Nucmer(reference, query)

        actual = nucmer.generate_command()
        expected = ["nucmer", "--prefix", "out", str(reference), str(query)]

        assert actual == expected

    def test_generateCommand_extraParamsIncluded(self):
        reference = REF
        query = QUERY
        extra_params = "--maxmatch --forward"
        nucmer = Nucmer(reference, query, extra_params=extra_params)

        actual = nucmer.generate_command()
        expected = [
            "nucmer",
            extra_params,
            "--prefix",
            "out",
            str(reference),
            str(query),
        ]

        assert actual == expected

    def test_run_twoTestCases_returnsOkExitCode(self):
        reference = REF
        query = QUERY
        prefix = f"/tmp/{str(uuid.uuid4())}"
        extra_params = "--maxmatch"
        nucmer = Nucmer(reference, query, prefix=prefix, extra_params=extra_params)
        result = nucmer.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        Path(prefix + ".delta").unlink()


class TestDeltaFilter:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        deltafile = Path("foo")

        with pytest.raises(DeltaFilterError):
            deltafilter = DeltaFilter(deltafile)

    def test_generateCommand_defaultArgsHasNoExtraParamsAndHasOutAsPrefix(self):
        deltafile = DELTA
        deltafilter = DeltaFilter(deltafile)

        actual = deltafilter.generate_command()
        expected = ["delta-filter", str(deltafile)]

        assert actual == expected

    def test_generateCommand_extraParamsIncluded(self):
        deltafile = DELTA
        extra_params = "-1"
        deltafilter = DeltaFilter(deltafile, extra_params=extra_params)

        actual = deltafilter.generate_command()
        expected = ["delta-filter", extra_params, str(deltafile)]

        assert actual == expected

    def test_run_realDeltaFileReturnsOkExitCode(self):
        deltafile = DELTA
        extra_params = "-1"
        deltafilter = DeltaFilter(deltafile, extra_params=extra_params)
        result = deltafilter.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        actual = result.stdout.decode()
        expected = """test_cases/ref.fa test_cases/query.fa
NUCMER
>ref query 85 84
1 85 1 84 2 2 0
39
0
"""

        assert actual == expected


class TestShowSnps:
    def test_init_withFakeFileRaisesFileNotFoundError(self):
        deltafile = Path("foo")

        with pytest.raises(ShowSnpsError):
            showsnps = ShowSnps(deltafile)

    def test_generateCommand_defaultArgsHasNoExtraParams(self):
        deltafile = DELTA1
        showsnps = ShowSnps(deltafile)

        actual = showsnps.generate_command()
        expected = ["show-snps", str(deltafile)]

        assert actual == expected

    def test_generateCommand_oppositeDefaultArgsAndAmbiguousMappingExtraParam(self):
        deltafile = DELTA1
        print_header = False
        indels = False
        context = 3
        extra_params = "-C"
        showsnps = ShowSnps(
            deltafile,
            context=context,
            print_header=print_header,
            indels=indels,
            extra_params=extra_params,
        )

        actual = showsnps.generate_command()
        expected = ["show-snps", "-C", f"-x {context}", "-H", "-I", str(deltafile)]

        assert actual == expected

    def test_run_realDeltaFileReturnsOkExitCodeAndExpectedOutput(self):
        deltafile = DELTA1
        context = 3
        extra_params = "-rlTC"
        showsnps = ShowSnps(deltafile, context=context, extra_params=extra_params)
        result = showsnps.run()

        actual = result.returncode
        expected = 0

        assert actual == expected

        actual = result.stdout.decode()
        expected = """test_cases/ref.fa test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
73\tT\tA\t72\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
"""

        assert actual == expected

    def test_toDataframe_emptyFileReturnsEmpty(self):
        snps = StringIO()

        actual = ShowSnps.to_dataframe(snps)

        assert actual.empty

    def test_toDataframe_invalidInputRaisesError(self):
        snps = StringIO("foo\nbar\nsome\nrandom\ntext\n")

        with pytest.raises(ValueError):
            ShowSnps.to_dataframe(snps)

    def test_toDataFrame_validInputReturnCorrectDataframe(self):
        snps = StringIO(
            """test_cases/ref.fa test_cases/query.fa
NUCMER

[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
73\tT\tA\t72\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
"""
        )

        actual = ShowSnps.to_dataframe(snps)
        expected = ShowSNPsDataframe(
            {
                "ref_pos": [39, 73],
                "ref_sub": ["G", "T"],
                "query_sub": [".", "A"],
                "query_pos": [38, 72],
                "nearest_mismatch": [34, 13],
                "nearest_end": [38, 13],
                "ref_len": [85, 85],
                "query_len": [84, 84],
                "ref_context": ["GTAGTAG", "GGATTGA"],
                "query_context": ["GTA.TAG", "GGAATGA"],
                "ref_strand": [1, 1],
                "query_strand": [1, 1],
                "ref_chrom": ["ref", "ref"],
                "query_chrom": ["query", "query"],
            }
        )

        assert actual.equals(expected)


class TestShowSNPsDataframe:
    ref_length = 20
    ref_pos = 3
    rev_ref_pos = 18
    query_length = 10
    query_pos = 8
    rev_query_pos = 3

    @staticmethod
    def create_showssnps_content(
            ref_strand_is_reversed: bool, query_strand_is_reversed: bool
    ) -> StringIO:
        ref_strand = -1 if ref_strand_is_reversed else 1
        query_strand = -1 if query_strand_is_reversed else 1
        return StringIO(
            f"whatever whatever\n"
            f"NUCMER\n"
            f"\n"
            f"[P1]	[SUB]	[SUB]	[P2]	[BUFF]	[DIST]	[LEN R]	[LEN Q]	[CTX R]	[CTX Q]	[FRM]	[TAGS]\n"
            f"{TestShowSNPsDataframe.ref_pos}	G	C	{TestShowSNPsDataframe.query_pos}	1	1	{TestShowSNPsDataframe.ref_length}"
            f"	{TestShowSNPsDataframe.query_length}	AA	AA	{ref_strand}	{query_strand}	ref	query"
        )

    @staticmethod
    def create_expected_showsnps_dataframe(
            ref_strand_is_reversed: bool, query_strand_is_reversed: bool
    ) -> ShowSNPsDataframe:
        if ref_strand_is_reversed:
            ref_pos = TestShowSNPsDataframe.rev_ref_pos
        else:
            ref_pos = TestShowSNPsDataframe.ref_pos

        if query_strand_is_reversed:
            query_pos = TestShowSNPsDataframe.rev_query_pos
        else:
            query_pos = TestShowSNPsDataframe.query_pos

        return ShowSNPsDataframe(
            {
                "ref_pos": [ref_pos],
                "ref_sub": ["G"],
                "query_sub": ["C"],
                "query_pos": [query_pos],
                "nearest_mismatch": [1],
                "nearest_end": [1],
                "ref_len": [TestShowSNPsDataframe.ref_length],
                "query_len": [TestShowSNPsDataframe.query_length],
                "ref_context": ["AA"],
                "query_context": ["AA"],
                "ref_strand": [1],
                "query_strand": [1],
                "ref_chrom": ["ref"],
                "query_chrom": ["query"],
            }
        )

    def test_translateToFWDStrand_queryStrandIsForwardRefStrandIsForwardReturnsPositionsWithNoChange(
            self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=False, query_strand_is_reversed=False
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=False, query_strand_is_reversed=False
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsReverseRefStrandIsForwardReturnsPositionsWithQueryChanged(
            self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=False, query_strand_is_reversed=True
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=False, query_strand_is_reversed=True
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsForwardRefStrandIsReverseReturnsPositionsWithRefChanged(
            self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=True, query_strand_is_reversed=False
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=True, query_strand_is_reversed=False
        )

        assert actual.equals(expected)

    def test_translateToFWDStrand_queryStrandIsReverseRefStrandIsReverseReturnsBothPositionsChanged(
            self
    ):
        showsnps_content = TestShowSNPsDataframe.create_showssnps_content(
            ref_strand_is_reversed=True, query_strand_is_reversed=True
        )
        df = ShowSnps.to_dataframe(showsnps_content)

        actual = df.translate_to_FWD_strand()
        expected = TestShowSNPsDataframe.create_expected_showsnps_dataframe(
            ref_strand_is_reversed=True, query_strand_is_reversed=True
        )

        assert actual.equals(expected)

    def test_makePosZeroBased_emptyReturnsEmpty(self):
        df = ShowSNPsDataframe()

        actual = df.make_pos_zero_based()
        expected = df

        assert actual.equals(expected)

    def test_makePosZeroBased_oneValueReturnsValueMinusOne(self):
        df = ShowSNPsDataframe({"ref_pos": [1], "query_pos": [8]})

        actual = df.make_pos_zero_based()
        expected = ShowSNPsDataframe({"ref_pos": [0], "query_pos": [7]})

        assert actual.equals(expected)

    def test_makePosZeroBased_twoValuesReturnsAllValuesMinusOne(self):
        df = ShowSNPsDataframe({"ref_pos": [1, 89], "query_pos": [8, 10]})

        actual = df.make_pos_zero_based()
        expected = ShowSNPsDataframe({"ref_pos": [0, 88], "query_pos": [7, 9]})

        assert actual.equals(expected)

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___absolute_path(self):
        filepath = "/home/leandro/git/pandora1_paper/analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = ShowSNPsDataframe._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query == "H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___relative_path(self):
        filepath = "analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = ShowSNPsDataframe._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query == "H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path(self):
        filepath = "CFT073_and_H131800734.snps_df.pickle"
        ref, query = ShowSNPsDataframe._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query == "H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path___trivial_names(self):
        filepath = "A_and_B.snps_df.pickle"
        ref, query = ShowSNPsDataframe._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "A" and query == "B"



class TestGetReportFromDeltaFile:
    def test_get_ref_and_query_aligned_bases_percentage(self):
        report = """
/home/leandro/git/pandora1_paper/snippy_debug/refs/063_STEC.ref_chrom.fa /home/leandro/git/pandora1_paper/snippy_debug/refs/CFT073.ref_chrom.fa
NUCMER

                               [REF]                [QRY]
[Sequences]
TotalSeqs                          1                    1
AlignedSeqs               1(100.00%)           1(100.00%)
UnalignedSeqs               0(0.00%)             0(0.00%)

[Bases]
TotalBases                   4905296              5155066
AlignedBases         4236038(86.36%)      4202046(81.51%)
UnalignedBases        669258(13.64%)       953020(18.49%)

[Alignments]
1-to-1                           228                  228
TotalLength                  4163094              4163231
AvgLength                   18259.18             18259.79
AvgIdentity                    98.94                98.94

M-to-M                           343                  343
TotalLength                  4310592              4310855
AvgLength                   12567.32             12568.09
AvgIdentity                    98.83                98.83

[Feature Estimates]
Breakpoints                      684                  685
Relocations                       61                   51
Translocations                     0                    0
Inversions                        34                   26

Insertions                       252                  240
InsertionSum                  757221              1014066
InsertionAvg                 3004.85              4225.27

TandemIns                          2                    1
TandemInsSum                     171                   95
TandemInsAvg                   85.50                95.00

[SNPs]
TotalSNPs                      40615                40615
GC                       1048(2.58%)          1034(2.55%)
GT                       1570(3.87%)          1552(3.82%)
GA                      7361(18.12%)         7508(18.49%)
AC                       1545(3.80%)          1537(3.78%)
AG                      7508(18.49%)         7361(18.12%)
AT                       1365(3.36%)          1338(3.29%)
TC                      7295(17.96%)         7462(18.37%)
TA                       1338(3.29%)          1365(3.36%)
TG                       1552(3.82%)          1570(3.87%)
CA                       1537(3.78%)          1545(3.80%)
CG                       1034(2.55%)          1048(2.58%)
CT                      7462(18.37%)         7295(17.96%)

TotalGSNPs                     15112                15112
CG                        254(1.68%)           295(1.95%)
CA                        559(3.70%)           575(3.80%)
CT                      2800(18.53%)         2775(18.36%)
AG                      2924(19.35%)         2814(18.62%)
AT                        459(3.04%)           468(3.10%)
AC                        575(3.80%)           559(3.70%)
GA                      2814(18.62%)         2924(19.35%)
GT                        607(4.02%)           582(3.85%)
GC                        295(1.95%)           254(1.68%)
TC                      2775(18.36%)         2800(18.53%)
TG                        582(3.85%)           607(4.02%)
TA                        468(3.10%)           459(3.04%)

TotalIndels                     2157                 2157
G.                       256(11.87%)          265(12.29%)
A.                       273(12.66%)          320(14.84%)
T.                       243(11.27%)          280(12.98%)
C.                       246(11.40%)          274(12.70%)
.A                       320(14.84%)          273(12.66%)
.G                       265(12.29%)          256(11.87%)
.T                       280(12.98%)          243(11.27%)
.C                       274(12.70%)          246(11.40%)

TotalGIndels                     172                  172
C.                        21(12.21%)            17(9.88%)
A.                        26(15.12%)           34(19.77%)
G.                        25(14.53%)             6(3.49%)
T.                        23(13.37%)           20(11.63%)
.C                         17(9.88%)           21(12.21%)
.T                        20(11.63%)           23(13.37%)
.A                        34(19.77%)           26(15.12%)
.G                          6(3.49%)           25(14.53%)

"""
        actual = GetReportFromDeltaFile.get_ref_and_query_aligned_bases_percentage(report)
        expected = (86.36, 81.51)
        assert actual == expected
