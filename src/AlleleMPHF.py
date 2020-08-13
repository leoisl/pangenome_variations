from typing import List
import logging

from src.Allele import Allele
from src.MPHF import MPHF
from src.VarifierDataframe import VarifierDataframe


class AlleleMPHF(MPHF):
    """
    This class just aggregates helper functions to be used by the pipeline, not tested
    """
    def __init__(self):
        super().__init__()

    def _add_variants_from_VarifierDataframe_core(self, ref: str, query: str, snps_df: VarifierDataframe):
        for ref_allele, query_allele in Allele.get_alleles_from_VarifierDataframe(ref, query, snps_df):
            self.add_object(ref_allele)
            self.add_object(query_allele)

    def _add_variants_from_VarifierDataframe_filepath(self, VarifierDataframe_filepath):
        ref, query = VarifierDataframe.get_ref_and_query_from_VarifierDataframe_filepath(
            VarifierDataframe_filepath)
        snps_df = VarifierDataframe.load_pickled(VarifierDataframe_filepath)
        self._add_variants_from_VarifierDataframe_core(ref, query, snps_df)

    @staticmethod
    def build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths: List[str]) -> "AlleleMPHF":
        logging.info("Building AlleleMPHF from list of VarifierDataframes...")
        allele_mphf = AlleleMPHF()
        for snps_df_filepath in snps_dfs_filepaths:
            logging.info(f"Adding {snps_df_filepath}...")
            allele_mphf._add_variants_from_VarifierDataframe_filepath(snps_df_filepath)
        logging.info("Building AlleleMPHF from list of VarifierDataframes - Done")
        return allele_mphf
