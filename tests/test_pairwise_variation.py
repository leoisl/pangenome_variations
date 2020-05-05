from unittest import TestCase
from unittest.mock import Mock

from src.AlleleMPHF import AlleleMPHF
from src.PairwiseVariation import PairwiseVariation, AlleleMPHFNotSetException


class TestPairwiseVariation(TestCase):
    def test___constructor___ordered_alleles(self):
        allele_1 = 1
        allele_2 = 2
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1_id == allele_1
        assert pairwise_variation.allele_2_id == allele_2
        assert pairwise_variation.original_ref_allele_id == allele_1
        assert pairwise_variation.original_query_allele_id == allele_2

    def test___constructor___unordered_alleles(self):
        allele_1 = 2
        allele_2 = 1
        pairwise_variation = PairwiseVariation(allele_1, allele_2)
        assert pairwise_variation.allele_1_id == allele_2
        assert pairwise_variation.allele_2_id == allele_1
        assert pairwise_variation.original_ref_allele_id == allele_1
        assert pairwise_variation.original_query_allele_id == allele_2

    def test___equality___equalPairwiseVariations(self):
        allele_1 = 2
        allele_2 = 1
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert pairwise_variation_1 == pairwise_variation_2

    def test___equality___differentPairwiseVariations(self):
        allele_1 = 1
        allele_2 = 2
        allele_3 = 3
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_1, allele_3)
        assert pairwise_variation_1 != pairwise_variation_2

    def test_hash_equalPairwiseVariations_returnsEqualHashValues(self):
        allele_1 = 1
        allele_2 = 2
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert hash(pairwise_variation_1) == hash(pairwise_variation_2)

    def test_hash_differentPairwiseVariations_returnsDifferentHashValues(self):
        allele_1 = 1
        allele_2 = 2
        allele_3 = 3
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_1, allele_3)
        assert hash(pairwise_variation_1) != hash(pairwise_variation_2)

    def test___share_allele___no_sharing(self):
        allele_1 = 1
        allele_2 = 2
        allele_3 = 3
        allele_4 = 4
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_3, allele_4)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == False

    def test___share_allele___one_allele_shared(self):
        allele_1 = 1
        allele_2 = 2
        allele_3 = 3
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_3)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == True

    def test___share_allele___both_alleles_shared(self):
        allele_1 = 1
        allele_2 = 2
        pairwise_variation_1 = PairwiseVariation(allele_1, allele_2)
        pairwise_variation_2 = PairwiseVariation(allele_2, allele_1)
        assert pairwise_variation_1.share_allele(pairwise_variation_2) == True

    def test___lt___allele_1_lt_other_allele_1(self):
        allele_1_1 = 1
        allele_2_1 = 10
        allele_1_2 = 2
        allele_2_2 = 10
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertTrue(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_eq_other_allele_1___allele_2_lt_other_allele_2(self):
        allele_1_1 = 1
        allele_2_1 = 5
        allele_1_2 = 1
        allele_2_2 = 10
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertTrue(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_eq_other_allele_1___allele_2_eq_other_allele_2(self):
        allele_1_1 = 1
        allele_2_1 = 2
        allele_1_2 = 1
        allele_2_2 = 2
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertFalse(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_eq_other_allele_1___allele_2_gt_other_allele_2(self):
        allele_1_1 = 1
        allele_2_1 = 10
        allele_1_2 = 1
        allele_2_2 = 5
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertFalse(pairwise_variation_1 < pairwise_variation_2)

    def test___lt___allele_1_gt_other_allele_1(self):
        allele_1_1 = 3
        allele_2_1 = 5
        allele_1_2 = 1
        allele_2_2 = 10
        pairwise_variation_1 = PairwiseVariation(allele_1_1, allele_2_1)
        pairwise_variation_2 = PairwiseVariation(allele_1_2, allele_2_2)
        self.assertFalse(pairwise_variation_1 < pairwise_variation_2)

    def test___allele_mphf_usage(self):
        allele_0_mock = Mock()
        allele_1_mock = Mock()
        allele_mphf = AlleleMPHF()
        allele_mphf.add_object(allele_0_mock)
        allele_mphf.add_object(allele_1_mock)
        pairwise_variation = PairwiseVariation(1, 0, allele_mphf)

        self.assertEqual(pairwise_variation.allele_1, allele_0_mock)
        self.assertEqual(pairwise_variation.allele_2, allele_1_mock)
        self.assertEqual(pairwise_variation.original_ref_allele, allele_1_mock)
        self.assertEqual(pairwise_variation.original_query_allele, allele_0_mock)

    def test___allele_mphf___not_set(self):
        pairwise_variation = PairwiseVariation(0, 0)
        with self.assertRaises(AlleleMPHFNotSetException):
            pairwise_variation.check_allele_mphf()