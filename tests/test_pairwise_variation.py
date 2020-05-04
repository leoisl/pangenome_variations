from unittest import TestCase

from src.Allele import Allele
from src.PairwiseVariation import PairwiseVariation


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
