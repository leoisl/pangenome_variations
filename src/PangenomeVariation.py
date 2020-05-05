from collections import defaultdict
from typing import Iterable, List

from src.Allele import Allele


class PangenomeVariation:
    """
    Class representing a pangenome variation or an equivalence class.
    Note that, differently from PairwiseVariation, this class needs to know its Alleles.
    """

    # Note: trivial method, not tested
    def __init__(self, id: int, alleles: Iterable[Allele]):
        self._id = id
        self._alleles = sorted(list(set(alleles)))
        self._unique_allele_sequences = self._get_unique_allele_sequences()

    @property
    def id(self) -> int:
        return self._id

    @property
    def alleles(self) -> List[Allele]:
        return self._alleles

    @property
    def unique_allele_sequences(self) -> List[str]:
        return self._unique_allele_sequences

    # Note: trivial method, not tested
    def __eq__(self, other: object):
        if isinstance(other, PangenomeVariation):
            return self.id == other.id and self.alleles == other.alleles
        else:
            return False

    def _get_unique_allele_sequences(self) -> List[str]:
        """
        Get the set of different allele sequences in this Pangenome Variation.
        If this is a SNP A -> C, then we have 2 different allele sequences (A and C)
        If this Pangenome Variations has all possible SNPs, then we would have 4 different allele sequences (ACGT)
        """
        set_of_unique_allele_sequences = {allele.sequence for allele in self.alleles}
        unique_allele_sequences = sorted(list(set_of_unique_allele_sequences))
        return unique_allele_sequences

    def is_consistent(self) -> bool:
        genomes_to_chrom_and_pos = defaultdict(set)
        for allele in self.alleles:
            genomes_to_chrom_and_pos[allele.genome].add((allele.chrom, allele.pos))
        genomes_to_have_consistent_alleles = {
            genome: len(chrom_and_pos) == 1 for genome, chrom_and_pos in genomes_to_chrom_and_pos.items()
        }
        all_genomes_have_consistent_alleles = all(genomes_to_have_consistent_alleles.values())
        return all_genomes_have_consistent_alleles

    # Note: trivial getters, not tested:
    def get_number_of_alleles(self) -> int:
        return len(self.alleles)

    def get_allele_index(self, allele: Allele) -> int:
        return self.alleles.index(allele)

    def get_number_of_different_allele_sequences(self) -> int:
        return len(self.unique_allele_sequences)

    def get_allele_sequence_index(self, allele: Allele) -> int:
        return self.unique_allele_sequences.index(allele.sequence)

    def __repr__(self):
        return str(vars(self))
