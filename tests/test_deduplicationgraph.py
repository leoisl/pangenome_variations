from unittest import TestCase
from unittest.mock import patch

from src.DeduplicationGraph import DeduplicationGraph, InconsistentPairwiseVariation


class TestDeduplicationGraph(TestCase):
    def test____index_pairwise_variation___single_variation(self, *mocks):
        deduplication_graph = DeduplicationGraph(5, [(0, 1)])
        actual_indexing = deduplication_graph.allele_to_pairwise_variations
        expected_indexing = [{0}, {0}, set(), set(), set()]
        self.assertListEqual(actual_indexing, expected_indexing)

    def test____index_pairwise_variation___single_variation___incorrectly_sorted(self, *mocks):
        with self.assertRaises(InconsistentPairwiseVariation):
            DeduplicationGraph(5, [(1, 0)])

    def test____index_pairwise_variation___single_variation___incorrectly_sorted___variation_with_same_allele(self, *mocks):
        with self.assertRaises(InconsistentPairwiseVariation):
            DeduplicationGraph(5, [(1, 1)])

    def test____index_pairwise_variation___single_variation___allele_does_not_exist(self, *mocks):
        with self.assertRaises(IndexError):
            DeduplicationGraph(5, [(1, 5)])

    def test____index_pairwise_variation___several_variations(self, *mocks):
        deduplication_graph = DeduplicationGraph(10, [
            # first a triangle
            (0, 1), (0, 2), (1, 2),
            # then two pairs connected by a node
            (3, 4), (5, 6), (3, 6)
        ])

        actual_indexing = deduplication_graph.allele_to_pairwise_variations
        expected_indexing = [{0, 1}, {0, 2}, {1, 2}, {3, 5}, {3}, {4}, {4, 5}, set(), set(), set()]

        self.assertListEqual(actual_indexing, expected_indexing)

    def test____index_pairwise_variation___several_variations___last_is_inconsistent(self, *mocks):
        with self.assertRaises(InconsistentPairwiseVariation):
            DeduplicationGraph(10, [
                # first a triangle
                (0, 1), (0, 2), (1, 2),
                # then two pairs connected by a node
                (3, 4), (5, 6), (3, 6),
                (1, 0)
            ])


    def test___get_pangenome_variations_defined_by_allele_ids___totally_unconnected_graph(self):
        deduplication_graph = DeduplicationGraph(6, [(0, 1), (2, 3), (4, 5)])
        actual_pangenome_variations = deduplication_graph.get_pangenome_variations_defined_by_allele_ids()
        expected_pangenome_variations = [{0, 1}, {2, 3}, {4, 5}]
        self.assertListEqual(actual_pangenome_variations, expected_pangenome_variations)


    def test___get_pangenome_variations_defined_by_allele_ids___totally_connected_graph(self):
        deduplication_graph = DeduplicationGraph(5, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3),
                                                     (2, 4), (3, 4)])
        actual_pangenome_variations = deduplication_graph.get_pangenome_variations_defined_by_allele_ids()
        expected_pangenome_variations = [{0, 1, 2, 3, 4}]
        self.assertListEqual(actual_pangenome_variations, expected_pangenome_variations)

    def test___get_pangenome_variations_defined_by_allele_ids___isolated_node_and_a_path(self):
        deduplication_graph = DeduplicationGraph(10, [(0, 1), (2, 9), (3, 9), (3, 8), (4, 8), (4, 7), (5, 7), (5, 6)])
        actual_pangenome_variations = deduplication_graph.get_pangenome_variations_defined_by_allele_ids()
        expected_pangenome_variations = [{0, 1}, {2, 3, 4, 5, 6, 7, 8, 9}]
        self.assertListEqual(actual_pangenome_variations, expected_pangenome_variations)

    def test___get_pangenome_variations_defined_by_allele_ids___two_and_three_sized_components(self):
        deduplication_graph = DeduplicationGraph(10, [(0, 2), (3, 4), (1, 4)])
        actual_pangenome_variations = deduplication_graph.get_pangenome_variations_defined_by_allele_ids()
        expected_pangenome_variations = [{0, 2}, {1, 3, 4}]
        self.assertListEqual(actual_pangenome_variations, expected_pangenome_variations)