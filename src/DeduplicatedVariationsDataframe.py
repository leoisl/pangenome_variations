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

        for prefix in probe_prefixes:
            header = ProbeHeader(
                sample=row[f"{prefix}_genome"],
                chrom=row[f"{prefix}_chrom"],
                pos=row[f"{prefix}_pos"],
                ref_length=len(row["ref_allele"]),
                interval=ProbeInterval.from_string(row[f"{prefix}_probe_interval"]),
                pangenome_variation_id=row["pangenome_variation_id"],
                number_of_alleles=row["number_of_alleles"],
                allele_id=row[f"{prefix}_allele_id"],
                number_of_different_allele_sequences=row["number_of_different_allele_sequences"],
                allele_sequence_id=row[f"{prefix}_allele_sequence_id"],
                nb_of_samples=row["nb_of_samples"],
            )
            full_sequence = row[f"{prefix}_probe"]
            probes.append(Probe(header=header, full_sequence=full_sequence))

        return tuple(probes)
