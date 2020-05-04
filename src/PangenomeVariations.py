from typing import List, Set

from src.AlleleMPHF import AlleleMPHF
from src.PangenomeVariation import PangenomeVariation


# Note: trivial class, not tested
class PangenomeVariations:
    """
    Stores a list of PangenomeVariation
    """

    def __init__(self):
        self._pangenome_variations = []

    @property
    def pangenome_variations(self) -> List[PangenomeVariation]:
        return self._pangenome_variations

    def append(self, pangenome_variation: PangenomeVariation):
        self.pangenome_variations.append(pangenome_variation)

    def __eq__(self, other: object):
        if isinstance(other, PangenomeVariations):
            return self.pangenome_variations == other.pangenome_variations
        else:
            return False

    def __repr__(self):
        return str(vars(self))

    @staticmethod
    def build_from_pangenome_variations_defined_by_allele_ids(
            pangenome_variations_defined_by_allele_ids: List[Set[int]],
            allele_mphf: AlleleMPHF) -> "PangenomeVariations":
        pangenome_variations = PangenomeVariations()
        for pangenome_variation_index, alleles_indexes in enumerate(pangenome_variations_defined_by_allele_ids):
            alleles = [allele_mphf.get_object(allele_index) for allele_index in alleles_indexes]
            pangenome_variation = PangenomeVariation(pangenome_variation_index, alleles)
            pangenome_variations.append(pangenome_variation)
        return pangenome_variations
