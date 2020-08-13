from unittest import TestCase
import pandas as pd
from io import StringIO

from src.Allele import NotASNP, Allele


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
        assert (allele_1 < allele_2) == False
        assert (allele_1 <= allele_2) == True
        assert (allele_1 > allele_2) == False
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

    def test___get_alleles_from_VarifierDataframe___empty_df___no_alleles_returned(self):
        snps_df = pd.read_csv(StringIO(
            """ref_pos,ref_allele,query_allele,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
            """
        ))
        actual_alleles = list(Allele.get_alleles_from_VarifierDataframe("genome_1", "genome_2", snps_df))
        expected_alleles = []
        self.assertListEqual(actual_alleles, expected_alleles)

    def test___get_alleles_from_VarifierDataframe___one_SNP_in_df___two_alleles_returned(self):
        snps_df = pd.read_csv(StringIO(
            """ref_pos,ref_allele,query_allele,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
            1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
            """
        ))
        actual_alleles = list(Allele.get_alleles_from_VarifierDataframe("genome_1", "genome_2", snps_df))
        expected_alleles = [(Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "C"))]
        self.assertListEqual(actual_alleles, expected_alleles)

    def test___get_alleles_from_VarifierDataframe___one_SNP_variation_four_non_SNP_variation_in_df___two_alleles_returned(self):
        snps_df = pd.read_csv(StringIO(
            """ref_pos,ref_allele,query_allele,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
            1,,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
            1,A,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
            1,A,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
            1,AA,CC,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
            1,AA,C,2,0,0,0,0,ACGT,ACGT,1,1,chrom_1,chrom_2
            """
        ), na_filter=False)
        actual_alleles = list(Allele.get_alleles_from_VarifierDataframe("genome_1", "genome_2", snps_df))
        expected_alleles = [(Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "C"))]
        self.assertListEqual(actual_alleles, expected_alleles)

    def test___get_alleles_from_VarifierDataframe___four_SNP_variation_four_non_SNP_variation_in_df___only_SNPs_alleles_returned(self):
        snps_df = pd.read_csv(StringIO(
            """ref_pos,ref_allele,query_allele,query_pos,nearest_mismatch,nearest_end,ref_len,query_len,ref_context,query_context,ref_strand,query_strand,ref_chrom,query_chrom
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
        actual_alleles = list(Allele.get_alleles_from_VarifierDataframe("genome_1", "genome_2", snps_df))
        expected_alleles = [(Allele("genome_1", "chrom_10", 10, "G"), Allele("genome_2", "chrom_20", 20, "T")),
                            (Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "C")),
                            (Allele("genome_1", "chrom_1", 1, "A"), Allele("genome_2", "chrom_2", 2, "C")),
                            (Allele("genome_1", "chrom_100", 100, "A"), Allele("genome_2", "chrom_200", 200, "T"))]
        self.assertListEqual(actual_alleles, expected_alleles)