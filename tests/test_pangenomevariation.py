from unittest import TestCase
from src.PangenomeVariation import PangenomeVariation
from src.Allele import Allele

class TestPangenomeVariation(TestCase):
    def test____get_unique_allele_sequences___one_allele(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A")])
        actual = pangenome_variation.unique_allele_sequences
        expected = ["A"]
        self.assertListEqual(actual, expected)

    def test____get_unique_allele_sequences___two_unique_alleles_with_several_alleles(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"), ])
        actual = pangenome_variation.unique_allele_sequences
        expected = ["A", "G"]
        self.assertListEqual(actual, expected)

    def test____get_unique_allele_sequences___three_unique_alleles_with_several_alleles(self, *mocks):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "G"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "C"), ])
        actual = pangenome_variation.unique_allele_sequences
        expected = ["A", "C", "G"]
        self.assertListEqual(actual, expected)

    def test___is_consistent___two_genomes___single_chrom_and_pos_for_each_genome(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A")])
        self.assertTrue(pangenome_variation.is_consistent())

    def test___is_consistent___two_genomes___different_chrom_for_genome_1(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_2", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())

    def test___is_consistent___two_genomes___different_chrom_for_genome_2(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_1", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())

    def test___is_consistent___two_genomes___different_pos_for_genome_1(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 2, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())

    def test___is_consistent___two_genomes___different_pos_for_genome_2(self):
        pangenome_variation = PangenomeVariation(0, [Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 1, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A"),
                                                     Allele("genome_1", "chrom_1", 1, "A"),
                                                     Allele("genome_2", "chrom_2", 2, "A")])
        self.assertFalse(pangenome_variation.is_consistent())
