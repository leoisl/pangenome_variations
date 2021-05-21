from unittest import TestCase
from unittest.mock import patch
import pandas as pd

from src.DeduplicatedVariationsDataframe import DeduplicatedVariationsDataframe
from src.probe import Probe, ProbeHeader, ProbeInterval


class TestDeduplicatedVariationsDataframe(TestCase):
    def test___get_probes___empty_dataframe_returns_empty_panel(self):
        df = DeduplicatedVariationsDataframe()
        actual = df.get_probes()
        expected = ("", "")
        assert actual == expected

    def test___get_probes___invalid_dataframe_raises_error(self):
        df = DeduplicatedVariationsDataframe(
            {"ref_pos": [39, 73], "ref_allele": ["G", "T"], "query_len": [84, 84]}
        )

        with self.assertRaises(KeyError):
            df.get_probes()

    @patch.object(DeduplicatedVariationsDataframe, DeduplicatedVariationsDataframe._get_ref_and_query_probe.__name__,
                  side_effect=[("a", "b"), ("c", "d"), ("e", "f")])
    def test___get_probes(self, *mocks):
        df = DeduplicatedVariationsDataframe({"dummy": [0, 1, 2]})
        actual = df.get_probes()
        expected = ("a\nc\ne", "b\nd\nf")
        self.assertEqual(actual, expected)

    def test____get_ref_and_query_probe___SNP(self):
        row = pd.Series(data={
            "ref_genome": "ref_sample",
            "query_genome": "query_sample",
            "ref_allele": "T",
            "query_allele": "A",
            "ref_probe": "GGATTGA",
            "query_probe": "GGAATGA",
            "ref_probe_interval": "[3,4)",
            "query_probe_interval": "[3,4)",
            "ref_chrom": "1",
            "query_chrom": "1",
            "ref_pos": 73,
            "ref_original_pos": 73,
            "ref_original_strand": 1,
            "ref_len": 5000,
            "query_pos": 72,
            "query_original_pos": 322,
            "query_original_strand": -1,
            "query_len": 9999,
            "pangenome_variation_id": 42,
            "number_of_alleles": 5,
            "ref_allele_id": 2,
            "query_allele_id": 4,
            "number_of_different_allele_sequences": 10,
            "ref_allele_sequence_id": 5,
            "query_allele_sequence_id": 8,
            "nb_of_samples": 7,
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
                    allele_sequence_id=5,
                    nb_of_samples=7,
                    original_pos=73,
                    original_strand=1,
                    contig_length=5000
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
                    allele_sequence_id=8,
                    nb_of_samples=7,
                    original_pos=322,
                    original_strand=-1,
                    contig_length=9999
                ),
                full_sequence="GGAATGA"
            )
        )

        self.assertEqual(actual, expected)







    # these are indels, no need to be tested, but let's keep them here as they are possibly needed later
    # def test____get_ref_and_query_probe___indel(self):
    #     row = pd.Series(data={
    #         "ref_genome": "ref_sample",
    #         "query_genome": "query_sample",
    #         "ref_allele": "G",
    #         "query_allele": ".",
    #         "ref_context": "GTAGTAG",
    #         "query_context": "GTA.TAG",
    #         "ref_chrom": "1",
    #         "query_chrom": "1",
    #         "ref_pos": 39,
    #         "query_pos": 38,
    #         "pangenome_variation_id": 42,
    #         "number_of_alleles": 5,
    #         "ref_allele_id": 2,
    #         "query_allele_id": 4,
    #         "number_of_different_allele_sequences": 10,
    #         "ref_allele_sequence_id": 5,
    #         "query_allele_sequence_id": 8,
    #         "nb_of_samples": 7,
    #     })
    #     actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)
    #
    #     expected = (
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="ref_sample",
    #                 chrom="1",
    #                 pos=39,
    #                 ref_length=1,
    #                 interval=ProbeInterval(3, 4),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=2,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=5,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="GTAGTAG"
    #         ),
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="query_sample",
    #                 chrom="1",
    #                 pos=38,
    #                 ref_length=0,
    #                 interval=ProbeInterval(3, 3),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=4,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=8,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="GTATAG"
    #         )
    #     )
    #
    #     self.assertEqual(actual, expected)
    #
    # def test____get_ref_and_query_probe___indel_inverted(self):
    #     row = pd.Series(data={
    #         "ref_genome": "ref_sample",
    #         "query_genome": "query_sample",
    #         "query_allele": "G",
    #         "ref_allele": ".",
    #         "query_context": "GTAGTAG",
    #         "ref_context": "GTA.TAG",
    #         "ref_chrom": "1",
    #         "query_chrom": "1",
    #         "ref_pos": 39,
    #         "query_pos": 38,
    #         "pangenome_variation_id": 42,
    #         "number_of_alleles": 5,
    #         "ref_allele_id": 2,
    #         "query_allele_id": 4,
    #         "number_of_different_allele_sequences": 10,
    #         "ref_allele_sequence_id": 5,
    #         "query_allele_sequence_id": 8,
    #         "nb_of_samples": 7,
    #     })
    #     actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)
    #
    #     expected = (
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="ref_sample",
    #                 chrom="1",
    #                 pos=39,
    #                 ref_length=0,
    #                 interval=ProbeInterval(3, 3),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=2,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=5,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="GTATAG"
    #         ),
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="query_sample",
    #                 chrom="1",
    #                 pos=38,
    #                 ref_length=1,
    #                 interval=ProbeInterval(3, 4),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=4,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=8,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="GTAGTAG"
    #         )
    #     )
    #
    #     self.assertEqual(actual, expected)
    #
    # def test____get_ref_and_query_probe___indel_near_gene_start_truncated_left_flank(self):
    #     row = pd.Series(data={
    #         "ref_genome": "ref_sample",
    #         "query_genome": "query_sample",
    #         "ref_allele": "G",
    #         "query_allele": ".",
    #         "ref_context": "--AGTAG",
    #         "query_context": "-TA.TAG",
    #         "ref_chrom": "1",
    #         "query_chrom": "1",
    #         "ref_pos": 1,
    #         "query_pos": 1,
    #         "pangenome_variation_id": 42,
    #         "number_of_alleles": 5,
    #         "ref_allele_id": 2,
    #         "query_allele_id": 4,
    #         "number_of_different_allele_sequences": 10,
    #         "ref_allele_sequence_id": 5,
    #         "query_allele_sequence_id": 8,
    #         "nb_of_samples": 7,
    #     })
    #     actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)
    #
    #     expected = (
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="ref_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=1,
    #                 interval=ProbeInterval(1, 2),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=2,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=5,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="AGTAG"
    #         ),
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="query_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=0,
    #                 interval=ProbeInterval(2, 2),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=4,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=8,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="TATAG"
    #         )
    #     )
    #
    #     self.assertEqual(actual, expected)

    # this is an indel, no need to test
    # def test____get_ref_and_query_probe___indel_near_gene_end_truncated_right_flank(self):
    #     row = pd.Series(data={
    #         "ref_genome": "ref_sample",
    #         "query_genome": "query_sample",
    #         "ref_allele": "G",
    #         "query_allele": ".",
    #         "ref_context": "AAAGTA-",
    #         "query_context": "ATA.T--",
    #         "ref_chrom": "1",
    #         "query_chrom": "1",
    #         "ref_pos": 1,
    #         "query_pos": 1,
    #         "pangenome_variation_id": 42,
    #         "number_of_alleles": 5,
    #         "ref_allele_id": 2,
    #         "query_allele_id": 4,
    #         "number_of_different_allele_sequences": 10,
    #         "ref_allele_sequence_id": 5,
    #         "query_allele_sequence_id": 8,
    #         "nb_of_samples": 7,
    #     })
    #     actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)
    #
    #     expected = (
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="ref_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=1,
    #                 interval=ProbeInterval(3, 4),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=2,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=5,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="AAAGTA"
    #         ),
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="query_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=0,
    #                 interval=ProbeInterval(3, 3),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=4,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=8,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="ATAT"
    #         )
    #     )
    #
    #     self.assertEqual(actual, expected)

    # this is an indel, no need to test
    # def test____get_ref_and_query_probe___indel_at_gene_start_truncated_left_flank(self):
    #     row = pd.Series(data={
    #         "ref_genome": "ref_sample",
    #         "query_genome": "query_sample",
    #         "ref_allele": "G",
    #         "query_allele": ".",
    #         "ref_context": "---GTAG",
    #         "query_context": "---.TAG",
    #         "ref_chrom": "1",
    #         "query_chrom": "1",
    #         "ref_pos": 1,
    #         "query_pos": 1,
    #         "pangenome_variation_id": 42,
    #         "number_of_alleles": 5,
    #         "ref_allele_id": 2,
    #         "query_allele_id": 4,
    #         "number_of_different_allele_sequences": 10,
    #         "ref_allele_sequence_id": 5,
    #         "query_allele_sequence_id": 8,
    #         "nb_of_samples": 7,
    #     })
    #     actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)
    #
    #     expected = (
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="ref_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=1,
    #                 interval=ProbeInterval(0, 1),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=2,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=5,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="GTAG"
    #         ),
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="query_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=0,
    #                 interval=ProbeInterval(0, 0),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=4,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=8,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="TAG"
    #         )
    #     )
    #
    #     self.assertEqual(actual, expected)

    # this is an indel, no need to test
    # def test____get_ref_and_query_probe___indel_at_gene_end_truncated_right_flank(self):
    #     row = pd.Series(data={
    #         "ref_genome": "ref_sample",
    #         "query_genome": "query_sample",
    #         "ref_allele": "G",
    #         "query_allele": ".",
    #         "ref_context": "AAAG---",
    #         "query_context": "AAA.---",
    #         "ref_chrom": "1",
    #         "query_chrom": "1",
    #         "ref_pos": 1,
    #         "query_pos": 1,
    #         "pangenome_variation_id": 42,
    #         "number_of_alleles": 5,
    #         "ref_allele_id": 2,
    #         "query_allele_id": 4,
    #         "number_of_different_allele_sequences": 10,
    #         "ref_allele_sequence_id": 5,
    #         "query_allele_sequence_id": 8,
    #         "nb_of_samples": 7,
    #     })
    #     actual = DeduplicatedVariationsDataframe._get_ref_and_query_probe(row)
    #
    #     expected = (
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="ref_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=1,
    #                 interval=ProbeInterval(3, 4),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=2,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=5,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="AAAG"
    #         ),
    #         Probe(
    #             header=ProbeHeader(
    #                 sample="query_sample",
    #                 chrom="1",
    #                 pos=1,
    #                 ref_length=0,
    #                 interval=ProbeInterval(3, 3),
    #                 pangenome_variation_id=42,
    #                 number_of_alleles=5,
    #                 allele_id=4,
    #                 number_of_different_allele_sequences=10,
    #                 allele_sequence_id=8,
    #                 nb_of_samples=7,
    #             ),
    #             full_sequence="AAA"
    #         )
    #     )
    #
    #     self.assertEqual(actual, expected)

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
