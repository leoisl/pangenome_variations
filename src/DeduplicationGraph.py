from typing import List, Tuple, Set, Generator

import networkx as nx


class InconsistentPairwiseVariation(Exception):
    pass


class DeduplicationGraph:
    def __init__(self, number_of_alleles: int, pairwise_variation_id_to_alleles_id: List[Tuple[int, int]]):
        self._number_of_alleles = number_of_alleles
        self._pairwise_variation_id_to_alleles_id = pairwise_variation_id_to_alleles_id
        self._graph = nx.Graph()
        self._build_graph_nodes()
        self._index()
        self._build_edges()

    @property
    def number_of_alleles(self) -> int:
        return self._number_of_alleles

    @property
    def pairwise_variation_id_to_alleles_id(self) -> List[Tuple[int, int]]:
        return self._pairwise_variation_id_to_alleles_id

    @property
    def number_of_pairwise_variations(self) -> int:
        return len(self.pairwise_variation_id_to_alleles_id)

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    @property
    def nodes(self):
        return self.graph.nodes

    @property
    def edges(self):
        return self.graph.edges

    @property
    def allele_to_pairwise_variations(self) -> List[Set[int]]:
        return self._allele_to_pairwise_variations

    def _build_graph_nodes(self):
        self._graph.add_nodes_from(range(self.number_of_pairwise_variations))

    def _index(self):
        self._allele_to_pairwise_variations = [set() for _ in range(self.number_of_alleles)]
        for pairwise_variation_id, (allele_id_1, allele_id_2) in enumerate(self.pairwise_variation_id_to_alleles_id):
            if allele_id_1 >= allele_id_2:
                raise InconsistentPairwiseVariation(
                    f"Pairwise variation id {pairwise_variation_id} not correctly sorted: ({allele_id_1}, {allele_id_2})")
            self._allele_to_pairwise_variations[allele_id_1].add(pairwise_variation_id)
            self._allele_to_pairwise_variations[allele_id_2].add(pairwise_variation_id)

    def _add_edge(self, variant_1, variant_2) -> None:
        self._graph.add_edge(variant_1, variant_2)

    def _build_edges(self) -> None:
        """
        Note: Should be called after all nodes are added.
        """
        for pairwise_variations in self.allele_to_pairwise_variations:
            if len(pairwise_variations) > 1:
                # connect the variations with a path
                pairwise_variations_as_list = list(pairwise_variations)
                for pairwise_variation_1, pairwise_variation_2 in \
                        zip(pairwise_variations_as_list, pairwise_variations_as_list[1:]):
                    self._add_edge(pairwise_variation_1, pairwise_variation_2)

    def _get_connected_components(self) -> Generator[Set[int], None, None]:
        return nx.connected_components(self.graph)

    def get_pangenome_variations_defined_by_allele_ids(self) -> List[Set[int]]:
        pangenome_variations_defined_by_allele_ids = []

        connected_components = self._get_connected_components()
        for connected_component_index, connected_component in enumerate(connected_components):
            allele_ids_in_connected_component = []
            for pairwise_variation_id in connected_component:
                allele_ids_in_connected_component.extend(
                    self.pairwise_variation_id_to_alleles_id[pairwise_variation_id])
            pangenome_variation = set(allele_ids_in_connected_component)
            pangenome_variations_defined_by_allele_ids.append(pangenome_variation)

        return pangenome_variations_defined_by_allele_ids

    def __repr__(self):
        return str(vars(self))
