from typing import Tuple

import pandas as pd

from src.probe import Probe, ProbeHeader, ProbeInterval


class DeduplicatedVariationsDataframe(pd.DataFrame):
    def get_probes(self) -> Tuple[str, str]:
        ref_probes = []
        query_probes = []

        for index, row in self.iterrows():
            ref_probe, query_probe = self._get_ref_and_query_probe(row)
            ref_probes.append(str(ref_probe))
            query_probes.append(str(query_probe))

        return (
            "\n".join(probe for probe in ref_probes if probe),
            "\n".join(probe for probe in query_probes if probe),
        )

    @property
    def _constructor(self):
        return DeduplicatedVariationsDataframe

    @staticmethod
    def _get_ref_and_query_probe(row: pd.Series) -> Tuple[Probe, ...]:
        probes = []
        probe_prefixes = ["ref", "query"]

        # TODO: is this really " - 1" or len(ref/query_sub)??
        flank_width = int((len(row[f"{probe_prefixes[0]}_context"]) - 1) / 2)

        for prefix in probe_prefixes:
            core_sequence = row[f"{prefix}_sub"].replace(".", "")
            left_flank = row[f"{prefix}_context"][:flank_width].replace("-", "")
            right_flank = row[f"{prefix}_context"][flank_width + 1:].replace("-", "")
            call_start_idx = len(left_flank)
            call_end_idx = call_start_idx + len(core_sequence)
            header = ProbeHeader(
                sample=row[f"{prefix}_genome"],
                chrom=row[f"{prefix}_chrom"],
                pos=row[f"{prefix}_pos"],
                ref_length=len(core_sequence),
                interval=ProbeInterval(call_start_idx, call_end_idx),
                pangenome_variation_id=row["pangenome_variation_id"],
                number_of_alleles=row["number_of_alleles"],
                allele_id=row[f"{prefix}_allele_id"],
                number_of_different_allele_sequences=row["number_of_different_allele_sequences"],
                allele_sequence_id=row[f"{prefix}_allele_sequence_id"],
                nb_of_samples=row["nb_of_samples"],
            )
            full_sequence = left_flank + core_sequence + right_flank
            probes.append(Probe(header=header, full_sequence=full_sequence))

        return tuple(probes)
