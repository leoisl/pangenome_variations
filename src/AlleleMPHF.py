from typing import List

from src.Allele import Allele
from src.MPHF import MPHF
from src.Utils import Utils
from src.mummer import ShowSNPsDataframe


class AlleleMPHF(MPHF):
    def __init__(self):
        super().__init__()

    # helpers
    def _add_variants_from_ShowSNPsDataframe_core(self, ref: str, query: str, snps_df: ShowSNPsDataframe):
        for ref_allele, query_allele in Allele.get_alleles_from_ShowSNPsDataframe(ref, query, snps_df):
            self.add_object(ref_allele)
            self.add_object(query_allele)

    def _add_variants_from_ShowSNPsDataframe_filepath(self, ShowSNPsDataframe_filepath):
        ref, query = Utils._get_ref_and_query_from_ShowSNPsDataframe_filepath(
            ShowSNPsDataframe_filepath)
        snps_df = Utils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        self._add_variants_from_ShowSNPsDataframe_core(ref, query, snps_df)

    @staticmethod
    def build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths: List[str]) -> "AlleleMPHF":
        allele_mphf = AlleleMPHF()
        for snps_df_filepath in snps_dfs_filepaths:
            allele_mphf._add_variants_from_ShowSNPsDataframe_filepath(snps_df_filepath)
        return allele_mphf
