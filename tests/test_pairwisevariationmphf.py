from unittest import TestCase

from src.PairwiseVariationMPHF import PairwiseVariationMPHF
from src.PairwiseVariation import PairwiseVariation

class TestPairwiseVariationMPHF(TestCase):
    def test___get_pairwise_variation_id_to_alleles_id(self):
        mphf = PairwiseVariationMPHF()
        pairwise_variation_1 = PairwiseVariation(0, 1)
        pairwise_variation_2 = PairwiseVariation(1, 2)
        pairwise_variation_3 = PairwiseVariation(2, 3)

        mphf.add_object(pairwise_variation_1)
        mphf.add_object(pairwise_variation_1)
        mphf.add_object(pairwise_variation_3)
        mphf.add_object(pairwise_variation_1)
        mphf.add_object(pairwise_variation_2)
        mphf.add_object(pairwise_variation_3)
        mphf.add_object(pairwise_variation_1)


        actual_pairwise_variation_id_to_alleles_id = mphf.get_pairwise_variation_id_to_alleles_id()
        expected_pairwise_variation_id_to_alleles_id = [(0,1), (2,3), (1,2)]

        self.assertListEqual(actual_pairwise_variation_id_to_alleles_id,
                             expected_pairwise_variation_id_to_alleles_id)