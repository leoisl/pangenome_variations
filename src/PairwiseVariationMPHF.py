import pickle
from typing import List, Tuple, BinaryIO, TextIO

from src.AlleleMPHF import AlleleMPHF
from src.MPHF import MPHF
from src.PairwiseVariation import PairwiseVariation
from src.Utils import Utils
from src.mummer import ShowSNPsDataframe


class PairwiseVariationMPHF(MPHF):
    def __init__(self):
        super().__init__()

    # helpers
    def _add_variants_from_ShowSNPsDataframe_core(self, ref: str, query: str, snps_df: ShowSNPsDataframe,
                                                  allele_mphf: AlleleMPHF):
        for pairwise_variation in PairwiseVariation.get_PairwiseVariation_from_ShowSNPsDataframe(ref, query, snps_df,
                                                                                                 allele_mphf):
            self.add_object(pairwise_variation)

    def _add_variants_from_ShowSNPsDataframe_filepath(self, ShowSNPsDataframe_filepath: str, allele_mphf: AlleleMPHF):
        ref, query = Utils._get_ref_and_query_from_ShowSNPsDataframe_filepath(
            ShowSNPsDataframe_filepath)
        snps_df = Utils._load_pickled_ShowSNPsDataframe(ShowSNPsDataframe_filepath)
        self._add_variants_from_ShowSNPsDataframe_core(ref, query, snps_df, allele_mphf)

    @staticmethod
    def build_from_list_of_snps_dfs_filepaths(snps_dfs_filepaths: List[str],
                                              allele_mphf_filepath: str) -> "PairwiseVariationMPHF":
        with open(allele_mphf_filepath, "rb") as allele_mphf_file:
            allele_mphf = AlleleMPHF.load(allele_mphf_file)

        pairwise_variation_mphf = PairwiseVariationMPHF()
        for snps_df_filepath in snps_dfs_filepaths:
            pairwise_variation_mphf._add_variants_from_ShowSNPsDataframe_filepath(snps_df_filepath, allele_mphf)
        return pairwise_variation_mphf

    def get_pairwise_variation_id_to_alleles_id(self) -> List[Tuple[int, int]]:
        pairwise_variation_id_to_alleles_id = []
        for pairwise_variation in self.id_to_object:
            pairwise_variation_id_to_alleles_id.append((pairwise_variation.allele_1_id, pairwise_variation.allele_2_id))
        return pairwise_variation_id_to_alleles_id

    def dump(self, file_with_nb_of_objects: TextIO, pickle_file: BinaryIO,
             pairwise_variation_id_to_alleles_id_file: BinaryIO):
        super().dump(file_with_nb_of_objects, pickle_file)
        pairwise_variation_id_to_alleles_id = self.get_pairwise_variation_id_to_alleles_id()
        pickle.dump(pairwise_variation_id_to_alleles_id, pairwise_variation_id_to_alleles_id_file)
