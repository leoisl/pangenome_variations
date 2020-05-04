from evaluate.deduplicate_pairwise_snps import DeduplicatePairwiseSNPsUtils, NotASNP, Allele, PairwiseVariation, \
    DeduplicationGraph, PangenomeVariation, PangenomeVariations, ConsistentPangenomeVariations, InconsistentPangenomeVariations, \
    DeduplicatedVariationsDataframe
import pandas as pd
from io import StringIO
from unittest.mock import patch, PropertyMock, Mock
from unittest import TestCase
from collections import defaultdict
from evaluate.probe import Probe, ProbeInterval, ProbeHeader


class TestDeduplicatedVariationsDataframe(TestCase):
    def test___get_probes___empty_dataframe_returns_empty_panel(self):
        df = DeduplicatedVariationsDataframe()
        actual = df.get_probes()
        expected = ("", "")
        assert actual == expected


    def test___get_probes___invalid_dataframe_raises_error(self):
        df = DeduplicatedVariationsDataframe(
            {"ref_pos": [39, 73], "ref_sub": ["G", "T"], "query_len": [84, 84]}
        )

        with self.assertRaises(KeyError):
            df.get_probes()


    @patch.object(DeduplicatedVariationsDataframe, DeduplicatedVariationsDataframe._get_ref_and_query_probe.__name__,
                  side_effect=[("a", "b"), ("c", "d"), ("e", "f")])
    def test___get_probes(self, *mocks):
        df = DeduplicatedVariationsDataframe({"dummy": [0,1,2]})
        actual = df.get_probes()
        expected = ("a\nc\ne", "b\nd\nf")
        self.assertEqual(actual, expected)


    def test____get_ref_and_query_probe___nucmer_line_1(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_sub": "G",
            "query_sub": ".",
            "ref_context": "GTAGTAG",
            "query_context": "GTA.TAG",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 39,
            "query_pos": 38,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=39,
                    ref_length=1,
                    interval=ProbeInterval(3, 4),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="GTAGTAG"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=38,
                    ref_length=0,
                    interval=ProbeInterval(3, 3),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="GTATAG"
            )
        )

        self.assertEqual(actual, expected)


    def test____get_ref_and_query_probe___nucmer_line_1_inverted(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "query_sub": "G",
            "ref_sub": ".",
            "query_context": "GTAGTAG",
            "ref_context": "GTA.TAG",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 39,
            "query_pos": 38,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=39,
                    ref_length=0,
                    interval=ProbeInterval(3, 3),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="GTATAG"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=38,
                    ref_length=1,
                    interval=ProbeInterval(3, 4),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="GTAGTAG"
            )
        )

        self.assertEqual(actual, expected)

    def test____get_ref_and_query_probe___nucmer_line_2(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_sub": "T",
            "query_sub": "A",
            "ref_context": "GGATTGA",
            "query_context": "GGAATGA",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 73,
            "query_pos": 72,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=73,
                    ref_length=1,
                    interval=ProbeInterval(3, 4),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="GGATTGA"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=72,
                    ref_length=1,
                    interval=ProbeInterval(3, 4),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="GGAATGA"
            )
        )

        self.assertEqual(actual, expected)


    def test____get_ref_and_query_probe___nucmer_line_probe_near_gene_start_truncated_left_flank(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_sub": "G",
            "query_sub": ".",
            "ref_context": "--AGTAG",
            "query_context": "-TA.TAG",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 1,
            "query_pos": 1,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=1,
                    ref_length=1,
                    interval=ProbeInterval(1, 2),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="AGTAG"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=1,
                    ref_length=0,
                    interval=ProbeInterval(2, 2),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="TATAG"
            )
        )

        self.assertEqual(actual, expected)


    def test____get_ref_and_query_probe___nucmer_line_probe_near_gene_end_truncated_right_flank(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_sub": "G",
            "query_sub": ".",
            "ref_context": "AAAGTA-",
            "query_context": "ATA.T--",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 1,
            "query_pos": 1,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=1,
                    ref_length=1,
                    interval=ProbeInterval(3, 4),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="AAAGTA"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=1,
                    ref_length=0,
                    interval=ProbeInterval(3, 3),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="ATAT"
            )
        )

        self.assertEqual(actual, expected)


    def test____get_ref_and_query_probe___nucmer_line_probe_at_gene_start_truncated_left_flank(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_sub": "G",
            "query_sub": ".",
            "ref_context": "---GTAG",
            "query_context": "---.TAG",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 1,
            "query_pos": 1,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=1,
                    ref_length=1,
                    interval=ProbeInterval(0, 1),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="GTAG"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=1,
                    ref_length=0,
                    interval=ProbeInterval(0, 0),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="TAG"
            )
        )

        self.assertEqual(actual, expected)


    def test____get_ref_and_query_probe___nucmer_line_probe_at_gene_end_truncated_right_flank(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_sub": "G",
            "query_sub": ".",
            "ref_context": "AAAG---",
            "query_context": "AAA.---",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 1,
            "query_pos": 1,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
        })
        actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)

        expected = (
            Probe(
                header=ProbeHeader(
                    sample="ref_sample",
                    chrom="1",
                    pos=1,
                    ref_length=1,
                    interval=ProbeInterval(3, 4),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=2,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=5
                ),
                full_sequence="AAAG"
            ),
            Probe(
                header=ProbeHeader(
                    sample="query_sample",
                    chrom="1",
                    pos=1,
                    ref_length=0,
                    interval=ProbeInterval(3, 3),
                    pangenome_variation_id=42,
                    number_of_alleles=5,
                    allele_id=4,
                    number_of_different_allele_sequences=10,
                    allele_sequence_id=8
                ),
                full_sequence="AAA"
            )
        )

        self.assertEqual(actual, expected)


    # TODO: these tests were removed, as we are not merging stuff back anymore - Put them back in if merge is enabled
    # TODO: these tests were removed, as we are not merging stuff back anymore - Put them back in if merge is enabled
    # TODO: these tests were removed, as we are not merging stuff back anymore - Put them back in if merge is enabled
    # def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForIndelInQuery(
    #     self
    # ):
    #     df = ShowSnps.to_dataframe(
    #         StringIO(
    #             """ref query
    # NUCMER
    #
    # [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    # 39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    # 73\tT\t.\t72\t13\t13\t85\t84\tGGATTTG\tGGA.TGA\t1\t1\tref\tquery
    # 74\tT\t.\t72\t13\t13\t85\t84\tGATTTGA\tGGA.TGA\t1\t1\tref\tquery
    # 75\tT\t.\t72\t13\t13\t85\t84\tATTTGAA\tGGA.TGA\t1\t1\tref\tquery
    # 79\tT\tA\t78\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    # """
    #         )
    #     )
    #
    #     actual = df.get_probes()
    #     expected_ref = str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
    #             full_sequence="GTAGTAG",
    #         )
    #     )
    #     expected_ref += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=73, interval=ProbeInterval(3, 6)),
    #             full_sequence="GGATTTGAA",
    #         )
    #     )
    #     expected_ref += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=79, interval=ProbeInterval(3, 4)),
    #             full_sequence="GGATTGA",
    #         )
    #     )
    #     expected_query = str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
    #             full_sequence="GTATAG",
    #         )
    #     )
    #     expected_query += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=72, interval=ProbeInterval(3, 3)),
    #             full_sequence="GGATGA",
    #         )
    #     )
    #     expected_query += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=78, interval=ProbeInterval(3, 4)),
    #             full_sequence="GGAATGA",
    #         )
    #     )
    #     expected = (str(expected_ref), str(expected_query))
    #
    #     assert actual == expected
    #
    # def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForIndelInRef(self):
    #     df = ShowSnps.to_dataframe(
    #         StringIO(
    #             """ref query
    # NUCMER
    #
    # [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    # 39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    # 72\t.\tT\t73\t13\t13\t85\t84\tGGA.TGA\tGGATTTG\t1\t1\tref\tquery
    # 72\t.\tT\t74\t13\t13\t85\t84\tGGA.TGA\tGATTTGA\t1\t1\tref\tquery
    # 72\t.\tT\t75\t13\t13\t85\t84\tGGA.TGA\tATTTGAA\t1\t1\tref\tquery
    # 79\tT\tA\t78\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    # """
    #         )
    #     )
    #
    #     actual = df.get_probes()
    #     expected_ref = str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
    #             full_sequence="GTAGTAG",
    #         )
    #     )
    #     expected_ref += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=72, interval=ProbeInterval(3, 3)),
    #             full_sequence="GGATGA",
    #         )
    #     )
    #     expected_ref += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=79, interval=ProbeInterval(3, 4)),
    #             full_sequence="GGATTGA",
    #         )
    #     )
    #     expected_query = str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
    #             full_sequence="GTATAG",
    #         )
    #     )
    #     expected_query += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=73, interval=ProbeInterval(3, 6)),
    #             full_sequence="GGATTTGAA",
    #         )
    #     )
    #     expected_query += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=78, interval=ProbeInterval(3, 4)),
    #             full_sequence="GGAATGA",
    #         )
    #     )
    #     expected = (str(expected_ref), str(expected_query))
    #
    #     assert actual == expected
    #
    # def test_makeTruthPanelFromSnpsDataframe_mergeConsecutiveRecordsForMnpInRef(self):
    #     df = ShowSnps.to_dataframe(
    #         StringIO(
    #             """ref query
    # NUCMER
    #
    # [P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[LEN R]\t[LEN Q]\t[CTX R]\t[CTX Q]\t[FRM]\t[TAGS]
    # 39\tG\t.\t38\t34\t38\t85\t84\tGTAGTAG\tGTA.TAG\t1\t1\tref\tquery
    # 72\tA\tT\t73\t13\t13\t85\t84\tGGAAGCA\tGGATTTG\t1\t1\tref\tquery
    # 73\tG\tT\t74\t13\t13\t85\t84\tGAAGCAA\tGATTTGA\t1\t1\tref\tquery
    # 74\tC\tT\t75\t13\t13\t85\t84\tAAGCAAA\tATTTGAA\t1\t1\tref\tquery
    # 79\tT\tA\t78\t13\t13\t85\t84\tGGATTGA\tGGAATGA\t1\t1\tref\tquery
    # """
    #         )
    #     )
    #
    #     actual = df.get_probes()
    #     expected_ref = str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=39, interval=ProbeInterval(3, 4)),
    #             full_sequence="GTAGTAG",
    #         )
    #     )
    #     expected_ref += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=72, interval=ProbeInterval(3, 6)),
    #             full_sequence="GGAAGCAAA",
    #         )
    #     )
    #     expected_ref += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="ref", pos=79, interval=ProbeInterval(3, 4)),
    #             full_sequence="GGATTGA",
    #         )
    #     )
    #     expected_query = str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=38, interval=ProbeInterval(3, 3)),
    #             full_sequence="GTATAG",
    #         )
    #     )
    #     expected_query += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=73, interval=ProbeInterval(3, 6)),
    #             full_sequence="GGATTTGAA",
    #         )
    #     )
    #     expected_query += "\n" + str(
    #         Probe(
    #             ProbeHeader(chrom="query", pos=78, interval=ProbeInterval(3, 4)),
    #             full_sequence="GGAATGA",
    #         )
    #     )
    #     expected = (str(expected_ref), str(expected_query))
    #
    #     assert actual == expected


class TestDeduplicatePairwiseSNPsUtils(TestCase):
    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___absolute_path(self):
        filepath = "/home/leandro/git/pandora1_paper/analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___relative_path(self):
        filepath = "analysis_output/recall/snps_dfs/CFT073/CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path(self):
        filepath = "CFT073_and_H131800734.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "CFT073" and query=="H131800734"

    def test___get_ref_and_query_from_ShowSNPsDataframe_filepath___local_path___trivial_names(self):
        filepath = "A_and_B.snps_df.pickle"
        ref, query = DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath(filepath)
        assert ref == "A" and query=="B"



class TestAllele(TestCase):
    def test_constructor_isSNP_constructorOK(self):
        Allele("genome_1", "chrom_1", 10, "A")
        self.assertTrue(True)

    def test_constructor_isNotASNP_constructorRaisesNotASNP(self):
        with self.assertRaises(NotASNP):
            Allele("genome_1", "chrom_1", 10, "AC")

    def test_constructor_isDeletion_constructorRaisesNotASNP(self):
        with self.assertRaises(NotASNP):
            Allele("genome_1", "chrom_1", 10, "")

    def test_hash_equalAlleles_returnsEqualHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        self.assertEqual(hash(allele_1), hash(allele_2))

    def test_hash_differentGenomes_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_1", 10, "A")
        assert hash(allele_1) != hash(allele_2)

    def test_hash_differentChroms_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_2", 10, "A")
        assert hash(allele_1) != hash(allele_2)

    def test_hash_differentPos_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 11, "A")
        assert hash(allele_1) != hash(allele_2)

    def test_hash_differentSeqs_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "C")
        assert hash(allele_1) != hash(allele_2)

    def test___relational_operators___equalAlleles(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == True
        assert (allele_1 != allele_2) == False
        assert (allele_1 <  allele_2) == False
        assert (allele_1 <= allele_2) == True
        assert (allele_1 >  allele_2) == False
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentGenomes___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentGenomes___first_is_larger(self):
        allele_1 = Allele("genome_2", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentChroms___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_2", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentChroms___first_is_larger(self):
        allele_1 = Allele("genome_1", "chrom_2", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentPos___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 11, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentPos___first_is_larger(self):
        allele_1 = Allele("genome_1", "chrom_1", 11, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True

    def test___relational_operators___differentSeqs___first_is_smaller(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "C")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == True
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
        assert (allele_1 >= allele_2) == False

    def test___relational_operators___differentSeqs___first_is_larger(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "C")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        assert (allele_1 == allele_2) == False
        assert (allele_1 != allele_2) == True
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == False
        assert (allele_1 > allele_2) == True
        assert (allele_1 >= allele_2) == True


class TestPairwiseVariation(TestCase):
    def test___constructor___ordered_alleles(self):
        allele_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2 = Allele("genome_2", "chrom_2", 20, "A")
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1 == allele_1
        assert pairwise_variation.allele_2 == allele_2
        assert pairwise_variation.original_ref_allele == allele_1
        assert pairwise_variation.original_query_allele == allele_2

    def test___constructor___unordered_alleles(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1 == allele_2
        assert pairwise_variation.allele_2 == allele_1
        assert pairwise_variation.original_ref_allele == allele_1
        assert pairwise_variation.original_query_allele == allele_2

    def test___equality___equalPairwiseVariations(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert pairwise_variation_1 == pairwise_variation_2

    def test___equality___differentPairwiseVariations(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        allele_3 = Allele("genome_1", "chrom_1", 11, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_1, allele_3)
        assert pairwise_variation_1 != pairwise_variation_2

    def test_hash_equalPairwiseVariations_returnsEqualHashValues(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert hash(pairwise_variation_1) == hash(pairwise_variation_2)

    def test_hash_differentPairwiseVariations_returnsDifferentHashValues(self):
        allele_1 = Allele("genome_2", "chrom_2", 20, "A")
        allele_2 = Allele("genome_1", "chrom_1", 10, "A")
        allele_3 = Allele("genome_1", "chrom_1", 11, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_1, allele_3)
        assert hash(pairwise_variation_1) != hash(pairwise_variation_2)

    def test___share_allele___no_sharing(self):
        allele_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2 = Allele("genome_1", "chrom_1", 2, "A")
        allele_3 = Allele("genome_1", "chrom_1", 3, "A")
        allele_4 = Allele("genome_1", "chrom_1", 4, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_3, allele_4)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == False

    def test___share_allele___one_allele_shared(self):
        allele_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2 = Allele("genome_1", "chrom_1", 2, "A")
        allele_3 = Allele("genome_1", "chrom_1", 3, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_3)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == True

    def test___share_allele___both_alleles_shared(self):
        allele_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2 = Allele("genome_1", "chrom_1", 2, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == True


    def test___lt___allele_1_lt_other_allele_1(self):
        allele_1_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_1 = Allele("genome_1", "chrom_1", 2, "A")
        allele_1_2 = Allele("genome_1", "chrom_1", 3, "A")
        allele_2_2 = Allele("genome_1", "chrom_1", 4, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertTrue(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_eq_other_allele_1___allele_2_lt_other_allele_2(self):
        allele_1_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_1 = Allele("genome_1", "chrom_1", 2, "A")
        allele_1_2 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_2 = Allele("genome_1", "chrom_1", 4, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertTrue(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_eq_other_allele_1___allele_2_eq_other_allele_2(self):
        allele_1_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_1 = Allele("genome_1", "chrom_1", 2, "A")
        allele_1_2 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_2 = Allele("genome_1", "chrom_1", 2, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertFalse(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_eq_other_allele_1___allele_2_gt_other_allele_2(self):
        allele_1_1 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_1 = Allele("genome_1", "chrom_1", 4, "A")
        allele_1_2 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_2 = Allele("genome_1", "chrom_1", 2, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertFalse(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_gt_other_allele_1(self):
        allele_1_1 = Allele("genome_1", "chrom_1", 10, "A")
        allele_2_1 = Allele("genome_1", "chrom_1", 2, "A")
        allele_1_2 = Allele("genome_1", "chrom_1", 1, "A")
        allele_2_2 = Allele("genome_1", "chrom_1", 2, "A")
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertFalse(pairwise_variation_1 < pairwise_variation_2)


class TestPangenomeVariation(TestCase):
    def test____get_unique_allele_sequences___one_allele(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A")])
        actual = pangenome_variation.unique_allele_sequences
        expected=["A"]
        self.assertListEqual(actual, expected)

    def test____get_unique_allele_sequences___two_unique_alleles_with_several_alleles(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),])
        actual = pangenome_variation.unique_allele_sequences
        expected=["A", "G"]
        self.assertListEqual(actual, expected)

    def test____get_unique_allele_sequences___three_unique_alleles_with_several_alleles(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "C"),])
        actual = pangenome_variation.unique_allele_sequences
        expected=["A", "C", "G"]
        self.assertListEqual(actual, expected)

    def test___is_consistent___two_genomes___single_chrom_and_pos_for_each_genome(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertTrue(pangenome_variation.is_consistent())

    def test___is_consistent___two_genomes___different_chrom_for_genome_1(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_2", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
    def test___is_consistent___two_genomes___different_chrom_for_genome_2(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_1", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
    def test___is_consistent___two_genomes___different_pos_for_genome_1(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 2, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
    def test___is_consistent___two_genomes___different_pos_for_genome_2(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())


class TestConsistentPangenomeVariations(TestCase):
    def setUp(self) -> None:
        self.dummy_consistent_pangenome_variations = ConsistentPangenomeVariations(PangenomeVariations())

    def test___constructor(self):
        # setup
        consistent_pangenome_variations = []
        alleles_to_consistent_pangenome_variations = {}
        for i in range(3):
            consistent_pangenome_variation = Mock()
            consistent_pangenome_variation.is_consistent.return_value = True
            consistent_pangenome_variation.alleles = [f"consistent_pangenome_variation_{i}.alleles"]
            alleles_to_consistent_pangenome_variations[f"consistent_pangenome_variation_{i}.alleles"] = consistent_pangenome_variation
            consistent_pangenome_variations.append(consistent_pangenome_variation)

        inconsistent_pangenome_variations = []
        for _ in range(3):
            inconsistent_pangenome_variation = Mock()
            inconsistent_pangenome_variation.is_consistent.return_value = False
            inconsistent_pangenome_variations.append(inconsistent_pangenome_variation)

        list_of_pangenome_variations = [
            consistent_pangenome_variations[0],
            consistent_pangenome_variations[1],
            inconsistent_pangenome_variations[0],
            inconsistent_pangenome_variations[1],
            consistent_pangenome_variations[2],
            inconsistent_pangenome_variations[2]
        ]
        pangenome_variations = PangenomeVariations()
        pangenome_variations._pangenome_variations = list_of_pangenome_variations
        actual_consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations)

        self.assertListEqual(actual_consistent_pangenome_variations.consistent_pangenome_variations, consistent_pangenome_variations)
        self.assertDictEqual(actual_consistent_pangenome_variations.alleles_to_consistent_pangenome_variations,
                             alleles_to_consistent_pangenome_variations)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1", "allele_2": "CPV1"}))
    def test___get_consistent_pangenome_variation___both_alleles_present_and_in_same_CPV(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        actual = self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)
        expected="CPV1"
        self.assertEqual(actual, expected)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1", "allele_2": "CPV2"}))
    def test___get_consistent_pangenome_variation___both_alleles_present_but_in_different_CPVs(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1"}))
    def test___get_consistent_pangenome_variation___only_first_allele_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_2": "CPV1"}))
    def test___get_consistent_pangenome_variation___only_second_allele_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations", new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None))
    def test___get_consistent_pangenome_variation___no_alleles_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        actual = self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)
        expected = None
        self.assertEqual(actual, expected)


    @patch.object(PangenomeVariation, "get_number_of_alleles", side_effect=[2,4])
    @patch.object(PangenomeVariation, "get_allele_index", side_effect=[1,0,2,3])
    @patch.object(PangenomeVariation, "get_number_of_different_allele_sequences", side_effect=[2,3])
    @patch.object(PangenomeVariation, "get_allele_sequence_index", side_effect=[0,0,2,1])
    @patch.object(ConsistentPangenomeVariations, "get_consistent_pangenome_variation",
                  side_effect=[None,PangenomeVariation(0, []), PangenomeVariation(1, []), None])
    @patch.object(PairwiseVariation, PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe.__name__,
                  return_value=[Mock(), Mock(), Mock(), Mock()])
    def test____enrich_ShowSNPsDataframe(self, *mocks):
        snps_df = pd.read_csv(StringIO(
"""dummy
0
1
2
3
"""
        ))
        actual = self.dummy_consistent_pangenome_variations._get_DeduplicatedVariationsDataframe("ref", "query", snps_df)

        expected=DeduplicatedVariationsDataframe(pd.read_csv(StringIO(
"""dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
0,ref,query,False,-1,-1,-1,-1,-1,-1,-1
1,ref,query,True,0,2,1,0,2,0,0
2,ref,query,True,1,4,2,3,3,2,1
3,ref,query,False,-1,-1,-1,-1,-1,-1,-1
"""
        )))

        self.assertTrue(actual.equals(expected))


    @patch.object(DeduplicatePairwiseSNPsUtils, DeduplicatePairwiseSNPsUtils._get_ref_and_query_from_ShowSNPsDataframe_filepath.__name__,
                  return_value=(None, None))
    @patch.object(DeduplicatePairwiseSNPsUtils, DeduplicatePairwiseSNPsUtils._load_pickled_ShowSNPsDataframe.__name__)
    @patch.object(ConsistentPangenomeVariations, ConsistentPangenomeVariations._get_DeduplicatedVariationsDataframe.__name__)
    def test___load_and_process_ShowSNPsDataframe(self, _enrich_ShowSNPsDataframe_mock, *other_mocks):
        _enrich_ShowSNPsDataframe_mock.return_value = pd.read_csv(StringIO(
"""dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
0,ref,query,False,-1,-1,-1,-1,-1,-1,-1
1,ref,query,True,0,2,1,0,2,0,0
2,ref,query,True,1,4,2,3,3,2,1
3,ref,query,False,-1,-1,-1,-1,-1,-1,-1
"""
        ))

        actual = self.dummy_consistent_pangenome_variations.build_DeduplicatedVariationsDataframe_from_ShowSNPsDataframe("dummy")
        expected=DeduplicatedVariationsDataframe(pd.read_csv(StringIO(
"""dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
1,ref,query,True,0,2,1,0,2,0,0
2,ref,query,True,1,4,2,3,3,2,1
"""
        )))

        self.assertTrue(actual.equals(expected))



class TestDeduplicationGraph(TestCase):
    def setUp(self):
        self.alleles = [Allele(f"genome_{i}", "chrom_1", 1, "A") for i in range(10)]
        self.pairwise_mutations = [PairwiseVariation(self.alleles[i*2], self.alleles[i*2+1]) for i in range(5)]
        self.deduplication_graph = DeduplicationGraph()

    def test____index_pairwise_variation___single_variation(self, *mocks):
        self.deduplication_graph._index_pairwise_variation(self.pairwise_mutations[0])
        expected_indexing = {self.alleles[0]: {self.pairwise_mutations[0]},
                             self.alleles[1]: {self.pairwise_mutations[0]}}
        self.assertDictEqual(self.deduplication_graph.allele_to_pairwise_variations, expected_indexing)

    def test____index_pairwise_variation___several_variations(self, *mocks):

        pairwise_variations = [
        # first a triangle
        PairwiseVariation(self.alleles[0], self.alleles[1]),
        PairwiseVariation(self.alleles[0], self.alleles[2]),
        PairwiseVariation(self.alleles[1], self.alleles[0]), # same as the first one
        PairwiseVariation(self.alleles[1], self.alleles[2]),
        PairwiseVariation(self.alleles[2], self.alleles[0]), # same as the second one
        PairwiseVariation(self.alleles[2], self.alleles[1]), # same as the fourth one
        # then two pairs connected by a node
        PairwiseVariation(self.alleles[3], self.alleles[4]),
        PairwiseVariation(self.alleles[5], self.alleles[6]),
        PairwiseVariation(self.alleles[3], self.alleles[6])]

        for variation in pairwise_variations:
            self.deduplication_graph._index_pairwise_variation(variation)

        expected_indexing = {self.alleles[0]: {pairwise_variations[0], pairwise_variations[1], pairwise_variations[2], pairwise_variations[4]},
                             self.alleles[1]: {pairwise_variations[0], pairwise_variations[2], pairwise_variations[3], pairwise_variations[5]},
                             self.alleles[2]: {pairwise_variations[1], pairwise_variations[3], pairwise_variations[4], pairwise_variations[5]},
                             self.alleles[3]: {pairwise_variations[6], pairwise_variations[8]},
                             self.alleles[4]: {pairwise_variations[6]},
                             self.alleles[5]: {pairwise_variations[7]},
                             self.alleles[6]: {pairwise_variations[7], pairwise_variations[8]}}

        self.assertDictEqual(self.deduplication_graph.allele_to_pairwise_variations, expected_indexing)


    def test___add_variants_from_ShowSNPsDataframe_core___empty_df___no_variations_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        graph_is_empty = self.deduplication_graph.graph.number_of_nodes()==0
        assert graph_is_empty

    def test___add_variants_from_ShowSNPsDataframe_core___one_variation_in_df___one_variation_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==1 and self.deduplication_graph.graph.has_node(pairwise_variation)


    def test___add_variants_from_ShowSNPsDataframe_core___one_SNP_variation_four_non_SNP_variation_in_df___only_SNP_variation_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ), na_filter=False)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==1 and self.deduplication_graph.graph.has_node(pairwise_variation)



    def test___add_variants_from_ShowSNPsDataframe_core___three_SNP_variation_four_non_SNP_variation_in_df___only_SNPs_variations_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
10,G,T,20,0,0,0,0,ACGT,ACGT,1,1,chrom_10,chrom_20
1,,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,AA,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
100,A,T,200,0,0,0,0,ACGT,ACGT,1,1,chrom_100,chrom_200
"""
        ), na_filter=False)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)


        pairwise_variation_1 = PairwiseVariation(Allele("genome_1", "chrom_10", 10, "G"),
                                                 Allele("genome_2", "chrom_20", 20, "T"))
        pairwise_variation_2 = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_3 = PairwiseVariation(Allele("genome_1", "chrom_100", 100, "A"),
                                                 Allele("genome_2", "chrom_200", 200, "T"))
        assert self.deduplication_graph.graph.number_of_nodes() == 3 and \
               self.deduplication_graph.graph.has_node(pairwise_variation_1) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_2) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_3)


    def test___add_variants_from_ShowSNPsDataframe_core___one_same_variation_in_three_dfs___one_variation_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)

        pairwise_variation = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                               Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==1 and self.deduplication_graph.graph.has_node(pairwise_variation)

    def test___add_variants_from_ShowSNPsDataframe_core___one_variation_in_three_dfs___three_variations_added(self):
        snps_df = pd.read_csv(StringIO(
"""ref_pos,ref_sub,query_sub,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
"""
        ))
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_1", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_3", "genome_2", snps_df)
        self.deduplication_graph._add_variants_from_ShowSNPsDataframe_core("genome_4", "genome_2", snps_df)

        pairwise_variation_1 = PairwiseVariation(Allele("genome_1", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_2 = PairwiseVariation(Allele("genome_3", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        pairwise_variation_3 = PairwiseVariation(Allele("genome_4", "chrom_1", 1, "A"),
                                                 Allele("genome_2", "chrom_2", 2, "C"))
        assert self.deduplication_graph.graph.number_of_nodes()==3 and \
               self.deduplication_graph.graph.has_node(pairwise_variation_1) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_2) and \
               self.deduplication_graph.graph.has_node(pairwise_variation_3)

    @patch.object(DeduplicationGraph, "allele_to_pairwise_variations", new_callable=PropertyMock, return_value={})
    def test___build_edges___no_shared_alleles___no_edge_built(self, *mocks):
        self.deduplication_graph._build_edges()
        self.assertEqual(self.deduplication_graph.graph.number_of_edges(), 0)

    @patch.object(DeduplicationGraph, "allele_to_pairwise_variations", new_callable=PropertyMock)
    def test___build_edges___complex_case(self, allele_to_pairwise_variations_mock):
        allele_to_pairwise_variations_mock.return_value = {
        "allele_0": {self.pairwise_mutations[0], self.pairwise_mutations[2]}, # add edges both ways
        "allele_1": {self.pairwise_mutations[2], self.pairwise_mutations[0]}, # add edges both ways
        "allele_2": {self.pairwise_mutations[0], self.pairwise_mutations[1]}, # add one way
        "allele_3": {self.pairwise_mutations[4], self.pairwise_mutations[2]}, # add the other way
        "allele_4": {self.pairwise_mutations[3], self.pairwise_mutations[3]}, # self loop, not added
    }
        self.deduplication_graph._build_edges()
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[0], self.pairwise_mutations[2]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[2], self.pairwise_mutations[0]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[1], self.pairwise_mutations[0]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[0], self.pairwise_mutations[1]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[4], self.pairwise_mutations[2]))
        self.assertTrue(self.deduplication_graph.graph.has_edge(self.pairwise_mutations[2], self.pairwise_mutations[4]))
        self.assertEqual(self.deduplication_graph.graph.number_of_edges(), 3)


    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___totally_unconnected_graph(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0]},
            {self.pairwise_mutations[1]},
            {self.pairwise_mutations[2]},
            {self.pairwise_mutations[3]},
            {self.pairwise_mutations[4]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations_defined_by_allele_ids()

        pangenome_variations_expected = PangenomeVariations()
        for i in range(5):
            pangenome_variations_expected.append(PangenomeVariation(i, [
                self.pairwise_mutations[i].allele_1, self.pairwise_mutations[i].allele_2]))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)


    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___totally_connected_graph(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0], self.pairwise_mutations[1], self.pairwise_mutations[2],
             self.pairwise_mutations[3], self.pairwise_mutations[4]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations_defined_by_allele_ids()

        pangenome_variations_expected = PangenomeVariations()
        pangenome_variations_expected.append(PangenomeVariation(0, self.alleles))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)

    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___isolated_node_and_a_path(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0], self.pairwise_mutations[1], self.pairwise_mutations[3], self.pairwise_mutations[4]},
            {self.pairwise_mutations[2]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations_defined_by_allele_ids()

        pangenome_variations_expected = PangenomeVariations()
        pangenome_variations_expected.append(PangenomeVariation(0, [
            self.pairwise_mutations[0].allele_1, self.pairwise_mutations[0].allele_2,
            self.pairwise_mutations[1].allele_1, self.pairwise_mutations[1].allele_2,
            self.pairwise_mutations[3].allele_1, self.pairwise_mutations[3].allele_2,
            self.pairwise_mutations[4].allele_1, self.pairwise_mutations[4].allele_2,
        ]))
        pangenome_variations_expected.append(PangenomeVariation(1, [
            self.pairwise_mutations[2].allele_1, self.pairwise_mutations[2].allele_2,
        ]))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)

    @patch.object(DeduplicationGraph, DeduplicationGraph._get_connected_components.__name__)
    def test___get_pangenome_variations___two_and_three_sized_components(self, get_connected_components_mock):
        get_connected_components_mock.return_value = [
            {self.pairwise_mutations[0], self.pairwise_mutations[1], self.pairwise_mutations[3]},
            {self.pairwise_mutations[2], self.pairwise_mutations[4]}
        ]
        pangenome_variations_actual = self.deduplication_graph.get_pangenome_variations_defined_by_allele_ids()

        pangenome_variations_expected = PangenomeVariations()
        pangenome_variations_expected.append(PangenomeVariation(0, [
            self.pairwise_mutations[0].allele_1, self.pairwise_mutations[0].allele_2,
            self.pairwise_mutations[1].allele_1, self.pairwise_mutations[1].allele_2,
            self.pairwise_mutations[3].allele_1, self.pairwise_mutations[3].allele_2,
        ]))
        pangenome_variations_expected.append(PangenomeVariation(1, [
            self.pairwise_mutations[2].allele_1, self.pairwise_mutations[2].allele_2,
            self.pairwise_mutations[4].allele_1, self.pairwise_mutations[4].allele_2,
        ]))

        self.assertEqual(pangenome_variations_actual, pangenome_variations_expected)