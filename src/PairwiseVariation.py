import functools
from typing import Generator

from src.Allele import Allele
from src.AlleleMPHF import AlleleMPHF
from src.VarifierDataframe import VarifierDataframe


class AlleleMPHFNotSetException(Exception):
    pass


@functools.total_ordering
class PairwiseVariation:
    """
    Note: Pairwise variation does not know alleles, only allele IDs.
    It can know alleles if the allele_mphf is given, which translate allele IDs to alleles.
    """
    def __init__(self, ref_allele_id: int, query_allele_id: int,
                 allele_mphf: AlleleMPHF = None):
        """
        :param allele_mphf: translator of allele IDs to alleles, making a PairwiseVariation actually know its alleles
        """
        # this ordering is done to facilitate the usage of this class, for the equal and hash functions
        self._allele_1_id = min(ref_allele_id, query_allele_id)
        self._allele_2_id = max(ref_allele_id, query_allele_id)
        self._original_ref_allele_id = ref_allele_id
        self._original_query_allele_id = query_allele_id
        self._allele_mphf = allele_mphf

    @property
    def allele_1_id(self) -> int:
        return self._allele_1_id

    @property
    def allele_2_id(self) -> int:
        return self._allele_2_id

    @property
    def original_ref_allele_id(self) -> int:
        return self._original_ref_allele_id

    @property
    def original_query_allele_id(self) -> int:
        return self._original_query_allele_id

    @property
    def allele_mphf(self) -> AlleleMPHF:
        return self._allele_mphf

    def check_allele_mphf(self):
        if self.allele_mphf is None:
            raise AlleleMPHFNotSetException()

    @property
    def allele_1(self) -> Allele:
        self.check_allele_mphf()
        return self.allele_mphf.get_object(self.allele_1_id)

    @property
    def allele_2(self) -> Allele:
        self.check_allele_mphf()
        return self.allele_mphf.get_object(self.allele_2_id)

    @property
    def original_ref_allele(self) -> Allele:
        self.check_allele_mphf()
        return self.allele_mphf.get_object(self.original_ref_allele_id)

    @property
    def original_query_allele(self) -> Allele:
        self.check_allele_mphf()
        return self.allele_mphf.get_object(self.original_query_allele_id)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PairwiseVariation):
            return (self.allele_1_id, self.allele_2_id) == (other.allele_1_id, other.allele_2_id)
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.allele_1_id, self.allele_2_id))

    def __lt__(self, other: object) -> bool:
        if isinstance(other, PairwiseVariation):
            return (self.allele_1_id, self.allele_2_id) < (other.allele_1_id, other.allele_2_id)
        else:
            raise TypeError()

    def __repr__(self):
        return str(vars(self))

    def share_allele(self, other: "PairwiseVariation") -> bool:
        return self.allele_1_id == other.allele_1_id or self.allele_1_id == other.allele_2_id or \
               self.allele_2_id == other.allele_1_id or self.allele_2_id == other.allele_2_id

    @staticmethod
    def get_PairwiseVariation_from_VarifierDataframe(ref: str, query: str, snps_df: VarifierDataframe,
                                                     allele_mphf: AlleleMPHF) -> Generator[
        "PairwiseVariation", None, None]:
        for ref_allele, query_allele in Allele.get_alleles_from_VarifierDataframe(ref, query, snps_df):
            ref_allele_id = allele_mphf.get_id(ref_allele)
            query_allele_id = allele_mphf.get_id(query_allele)
            yield PairwiseVariation(ref_allele_id, query_allele_id, allele_mphf)
