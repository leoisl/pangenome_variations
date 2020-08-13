import pickle
from typing import List, Tuple
import logging

from src.AlleleMPHF import AlleleMPHF
from src.MPHF import MPHF
from src.PairwiseVariation import PairwiseVariation
from src.VarifierDataframe import VarifierDataframe


class PairwiseVariationMPHF(MPHF):
    """
    This class mostly aggregates helper functions to be used by the pipeline, some functions are thus not tested
    """
    def __init__(self):
        super().__init__()

    # Note: not tested
    def _add_variants_from_VarifierDataframe_core(self, ref: str, query: str, snps_df: VarifierDataframe,
                                                  allele_mphf: AlleleMPHF):
        for pairwise_variation in PairwiseVariation.get_PairwiseVariation_from_VarifierDataframe(ref, query, snps_df,
                                                                                                 allele_mphf):
            self.add_object(pairwise_variation)

    # Note: not tested
    def _add_variants_from_VarifierDataframe_filepath(self, VarifierDataframe_filepath: str, allele_mphf: AlleleMPHF):
        ref, query = VarifierDataframe.get_ref_and_query_from_VarifierDataframe_filepath(
            VarifierDataframe_filepath)
        snps_df = VarifierDataframe.load_pickled(VarifierDataframe_filepath)
        self._add_variants_from_VarifierDataframe_core(ref, query, snps_df, allele_mphf)

    # Note: not tested
    @staticmethod
    def build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths: List[str],
                                              allele_mphf_filepath: str) -> "PairwiseVariationMPHF":
        logging.info("Building PairwiseVariationMPHF from list of VarifierDataframes...")

        logging.info("Loading AlleleMPHF...")
        allele_mphf = AlleleMPHF.load(allele_mphf_filepath)

        pairwise_variation_mphf = PairwiseVariationMPHF()
        for snps_df_filepath in snps_dfs_filepaths:
            logging.info(f"Adding {snps_df_filepath}...")
            pairwise_variation_mphf._add_variants_from_VarifierDataframe_filepath(snps_df_filepath, allele_mphf)
        logging.info("Building PairwiseVariationMPHF from list of VarifierDataframes - Done")

        return pairwise_variation_mphf

    def get_pairwise_variation_id_to_alleles_id(self) -> List[Tuple[int, int]]:
        pairwise_variation_id_to_alleles_id = []
        for pairwise_variation in self.id_to_object:
            pairwise_variation_id_to_alleles_id.append((pairwise_variation.allele_1_id, pairwise_variation.allele_2_id))
        return pairwise_variation_id_to_alleles_id

    # Note: not tested
    def save(self, file_with_nb_of_objects_filepath: str, pickle_filepath: str,
             pairwise_variation_id_to_alleles_id_filepath: str):
        super().save(file_with_nb_of_objects_filepath, pickle_filepath)

        pairwise_variation_id_to_alleles_id = self.get_pairwise_variation_id_to_alleles_id()
        with open(pairwise_variation_id_to_alleles_id_filepath, "wb") as pairwise_variation_id_to_alleles_id_filehandler:
            pickle.dump(pairwise_variation_id_to_alleles_id, pairwise_variation_id_to_alleles_id_filehandler)
