import functools
from typing import Tuple, Generator

from src.mummer import ShowSNPsDataframe


class NotASNP(Exception):
    pass


@functools.total_ordering
class Allele:
    def __init__(self, genome: str, chrom: str, pos: int, sequence: str):
        self._genome = genome
        self._chrom = chrom
        self._pos = pos  # TODO: do we need to care about strand when looking at pos?
        self._sequence = sequence

        # sequence is just used to ensure this is a SNP, for now
        is_snp = len(sequence) == 1
        if not is_snp:
            raise NotASNP()

    @property
    def genome(self) -> str:
        return self._genome

    @property
    def chrom(self) -> str:
        return self._chrom

    @property
    def pos(self) -> int:
        return self._pos

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def data_tuple(self) -> Tuple[str, str, int, str]:
        """
        :return: a tuple with all the data in the allele
        """
        return self.genome, self.chrom, self.pos, self.sequence

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Allele):
            return self.data_tuple == other.data_tuple
        else:
            return False

    def __hash__(self) -> int:
        return hash(self.data_tuple)

    def __lt__(self, other: object) -> bool:
        if isinstance(other, Allele):
            return self.data_tuple < other.data_tuple
        else:
            raise TypeError()

    def __repr__(self):
        return str(vars(self))

    @staticmethod
    def get_alleles_from_ShowSNPsDataframe(ref: str, query: str, snps_df: ShowSNPsDataframe) -> Generator[
        Tuple["Allele", "Allele"], None, None]:
        for ref_chrom, ref_pos, ref_sub, query_chrom, query_pos, query_sub in \
                zip(snps_df["ref_chrom"], snps_df["ref_pos"], snps_df["ref_sub"],
                    snps_df["query_chrom"], snps_df["query_pos"], snps_df["query_sub"]):
            is_snp = len(ref_sub) == 1 and len(query_sub) == 1
            if not is_snp:
                continue  # we just deal with SNPs as of now

            ref_allele = Allele(ref, ref_chrom, ref_pos, ref_sub)
            query_allele = Allele(query, query_chrom, query_pos, query_sub)
            yield ref_allele, query_allele
