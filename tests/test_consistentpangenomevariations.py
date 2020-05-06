from unittest import TestCase
from unittest.mock import Mock, PropertyMock, patch
from collections import defaultdict
import pandas as pd
from io import StringIO

from src.ConsistentPangenomeVariations import ConsistentPangenomeVariations, InconsistentPangenomeVariations
from src.DeduplicatedVariationsDataframe import DeduplicatedVariationsDataframe
from src.PangenomeVariations import PangenomeVariations
from src.PangenomeVariation import PangenomeVariation
from src.PairwiseVariation import PairwiseVariation
from src.mummer import ShowSNPsDataframe
from src.AlleleMPHF import AlleleMPHF


class TestConsistentPangenomeVariations(TestCase):
    def setUp(self) -> None:
        self.dummy_consistent_pangenome_variations = ConsistentPangenomeVariations(PangenomeVariations(), filter_for_biallelic=False)

    def test___constructor___filter_inconsistent_variations_out(self):
        # setup
        consistent_pangenome_variations = []
        alleles_to_consistent_pangenome_variations = {}
        for i in range(3):
            consistent_pangenome_variation = Mock()
            consistent_pangenome_variation.is_consistent.return_value = True
            consistent_pangenome_variation.alleles = [f"consistent_pangenome_variation_{i}.alleles"]
            alleles_to_consistent_pangenome_variations[
                f"consistent_pangenome_variation_{i}.alleles"] = consistent_pangenome_variation
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
        actual_consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations, filter_for_biallelic=False)

        self.assertListEqual(actual_consistent_pangenome_variations.consistent_pangenome_variations,
                             consistent_pangenome_variations)
        self.assertDictEqual(actual_consistent_pangenome_variations.alleles_to_consistent_pangenome_variations,
                             alleles_to_consistent_pangenome_variations)
        self.assertEqual(actual_consistent_pangenome_variations.number_of_pangenome_variations, 6)
        self.assertEqual(actual_consistent_pangenome_variations.number_of_consistent_pangenome_variations, 3)

    def test___constructor___filter_inconsistent_variations_out___keep_only_biallelic_ones(self):
        # setup
        biallelic_consistent_pangenome_variations = []
        non_biallelic_consistent_pangenome_variations = []
        alleles_to_biallelic_consistent_pangenome_variations = {}
        number_of_alleles_in_each_consistent_variation = [3, 1, 2, 4, 2]
        for consistent_variation_index, number_of_alleles in enumerate(number_of_alleles_in_each_consistent_variation):
            consistent_pangenome_variation = Mock(get_number_of_different_allele_sequences=Mock(return_value=number_of_alleles))
            consistent_pangenome_variation.is_consistent.return_value = True
            consistent_pangenome_variation.alleles = [
                f"consistent_pangenome_variation_{consistent_variation_index}.allele_{allele_index}"
                for allele_index in range(number_of_alleles)]

            if number_of_alleles == 2:
                for allele in consistent_pangenome_variation.alleles:
                    alleles_to_biallelic_consistent_pangenome_variations[allele] = consistent_pangenome_variation
                biallelic_consistent_pangenome_variations.append(consistent_pangenome_variation)
            else:
                non_biallelic_consistent_pangenome_variations.append(consistent_pangenome_variation)

        inconsistent_pangenome_variations = []
        for _ in range(3):
            inconsistent_pangenome_variation = Mock()
            inconsistent_pangenome_variation.is_consistent.return_value = False
            inconsistent_pangenome_variations.append(inconsistent_pangenome_variation)

        list_of_pangenome_variations = [
            biallelic_consistent_pangenome_variations[0],
            inconsistent_pangenome_variations[0],
            non_biallelic_consistent_pangenome_variations[0],
            inconsistent_pangenome_variations[1],
            inconsistent_pangenome_variations[2],
            biallelic_consistent_pangenome_variations[1],
            non_biallelic_consistent_pangenome_variations[1],
            non_biallelic_consistent_pangenome_variations[2],
        ]
        pangenome_variations = PangenomeVariations()
        pangenome_variations._pangenome_variations = list_of_pangenome_variations
        actual_consistent_pangenome_variations = ConsistentPangenomeVariations(pangenome_variations, filter_for_biallelic=True)

        self.assertListEqual(actual_consistent_pangenome_variations.consistent_pangenome_variations,
                             biallelic_consistent_pangenome_variations)
        self.assertDictEqual(actual_consistent_pangenome_variations.alleles_to_consistent_pangenome_variations,
                             alleles_to_biallelic_consistent_pangenome_variations)
        self.assertEqual(actual_consistent_pangenome_variations.number_of_pangenome_variations, 8)
        self.assertEqual(actual_consistent_pangenome_variations.number_of_consistent_pangenome_variations, 5)
        self.assertEqual(actual_consistent_pangenome_variations.number_of_consistent_biallelic_pangenome_variations, 2)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations",
                  new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1", "allele_2": "CPV1"}))
    def test___get_consistent_pangenome_variation___both_alleles_present_and_in_same_CPV(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        actual = self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)
        expected = "CPV1"
        self.assertEqual(actual, expected)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations",
                  new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1", "allele_2": "CPV2"}))
    def test___get_consistent_pangenome_variation___both_alleles_present_but_in_different_CPVs(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations",
                  new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_1": "CPV1"}))
    def test___get_consistent_pangenome_variation___only_first_allele_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations",
                  new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None, {"allele_2": "CPV1"}))
    def test___get_consistent_pangenome_variation___only_second_allele_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        with self.assertRaises(InconsistentPangenomeVariations):
            self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)

    @patch.object(ConsistentPangenomeVariations, "alleles_to_consistent_pangenome_variations",
                  new_callable=PropertyMock,
                  return_value=defaultdict(lambda: None))
    def test___get_consistent_pangenome_variation___no_alleles_present(self, *mocks):
        pairwise_variation = Mock(allele_1="allele_1", allele_2="allele_2")
        actual = self.dummy_consistent_pangenome_variations.get_consistent_pangenome_variation(pairwise_variation)
        expected = None
        self.assertEqual(actual, expected)

    @patch.object(PangenomeVariation, "get_number_of_alleles", side_effect=[2, 4])
    @patch.object(PangenomeVariation, "get_allele_index", side_effect=[1, 0, 2, 3])
    @patch.object(PangenomeVariation, "get_number_of_different_allele_sequences", side_effect=[2, 3])
    @patch.object(PangenomeVariation, "get_allele_sequence_index", side_effect=[0, 0, 2, 1])
    @patch.object(ConsistentPangenomeVariations, "get_consistent_pangenome_variation",
                  side_effect=[None, PangenomeVariation(0, []), PangenomeVariation(1, []), None])
    @patch.object(PairwiseVariation, PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe.__name__,
                  return_value=[Mock(), Mock(), Mock(), Mock()])
    def test____get_DeduplicatedVariationsDataframe(self, *mocks):
        snps_df = pd.read_csv(StringIO(
            """dummy
            0
            1
            2
            3
            """
        ))
        dummy_allele_mphf = AlleleMPHF()
        actual = self.dummy_consistent_pangenome_variations._get_DeduplicatedVariationsDataframe("ref", "query",
                                                                                                 snps_df, dummy_allele_mphf)

        expected = DeduplicatedVariationsDataframe(pd.read_csv(StringIO(
            """dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
            0,ref,query,False,-1,-1,-1,-1,-1,-1,-1
            1,ref,query,True,0,2,1,0,2,0,0
            2,ref,query,True,1,4,2,3,3,2,1
            3,ref,query,False,-1,-1,-1,-1,-1,-1,-1
            """
        )))

        self.assertTrue(actual.equals(expected))

    @patch.object(ShowSNPsDataframe,
                  ShowSNPsDataframe.get_ref_and_query_from_ShowSNPsDataframe_filepath.__name__,
                  return_value=(None, None))
    @patch.object(ShowSNPsDataframe, ShowSNPsDataframe.load_pickled.__name__)
    @patch.object(ConsistentPangenomeVariations,
                  ConsistentPangenomeVariations._get_DeduplicatedVariationsDataframe.__name__)
    def test___build_DeduplicatedVariationsDataframe_from_ShowSNPsDataframe(self, _enrich_ShowSNPsDataframe_mock, *other_mocks):
        _enrich_ShowSNPsDataframe_mock.return_value = pd.read_csv(StringIO(
            """dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
            0,ref,query,False,-1,-1,-1,-1,-1,-1,-1
            1,ref,query,True,0,2,1,0,2,0,0
            2,ref,query,True,1,4,2,3,3,2,1
            3,ref,query,False,-1,-1,-1,-1,-1,-1,-1
            """
        ))
        dummy_allele_mphf = AlleleMPHF()

        actual = self.dummy_consistent_pangenome_variations.build_DeduplicatedVariationsDataframe_from_ShowSNPsDataframe(
            "dummy", dummy_allele_mphf)
        expected = DeduplicatedVariationsDataframe(pd.read_csv(StringIO(
            """dummy,ref_genome,query_genome,present_in_a_consistent_pangenome_variation,pangenome_variation_id,number_of_alleles,ref_allele_id,query_allele_id,number_of_different_allele_sequences,ref_allele_sequence_id,query_allele_sequence_id
            1,ref,query,True,0,2,1,0,2,0,0
            2,ref,query,True,1,4,2,3,3,2,1
            """
        )))

        self.assertTrue(actual.equals(expected))


