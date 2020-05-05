from unittest import TestCase

from src.Allele import Allele
from src.AlleleMPHF import AlleleMPHF
from src.PangenomeVariation import PangenomeVariation
from src.PangenomeVariations import PangenomeVariations


class TestPangenomeVariations(TestCase):
    def test___build_from_pangenome_variations_defined_by_allele_ids(self):
        pangenome_variations_defined_by_allele_ids = [{0, 2, 3}, {1, 4}, {5, 6}]
        allele_mphf = AlleleMPHF()
        alleles = [Allele(str(i), str(i), i, "A") for i in range(7)]

        for allele in alleles:
            allele_mphf.add_object(allele)

        actual_pangenome_variations = PangenomeVariations.build_from_pangenome_variations_defined_by_allele_ids(
            pangenome_variations_defined_by_allele_ids, allele_mphf
        )
        expected_pangenome_variations = PangenomeVariations()
        expected_pangenome_variations.append(PangenomeVariation(0, [alleles[0], alleles[2], alleles[3]]))
        expected_pangenome_variations.append(PangenomeVariation(1, [alleles[1], alleles[4]]))
        expected_pangenome_variations.append(PangenomeVariation(2, [alleles[5], alleles[6]]))

        self.assertEqual(actual_pangenome_variations, expected_pangenome_variations)