import re
from typing import NamedTuple, Type, Optional

DELIM = ";"


class RegexError(Exception):
    pass


class ProbeInterval(NamedTuple):
    start: int = -1
    end: int = -1

    def __bool__(self) -> bool:
        return not self.is_null()

    def __str__(self) -> str:
        if self.is_null():
            return ""

        return f"[{self.start},{self.end})"

    def __len__(self) -> int:
        return self.end - self.start

    def is_null(self) -> bool:
        return self.start == -1 and self.end == -1

    interval_regex = re.compile(r"\[(\d+),(\d+)\)")

    @staticmethod
    def from_string(string: str) -> Optional["ProbeInterval"]:
        try:
            match = ProbeInterval.interval_regex.search(string)
            return ProbeInterval(int(match.group(1)), int(match.group(2)))
        except (TypeError, AttributeError):
            return None


class ProbeHeader:
    def __init__(
            self,
            sample: str = None,
            chrom: str = None,
            pos: int = None,
            ref_length: int = None,
            interval: ProbeInterval = None,
            svtype: str = None,
            gt_conf: float = None,
            coverage: float = None,
            pangenome_variation_id: int = None,
            number_of_alleles: int = None,
            allele_id: int = None,
            number_of_different_allele_sequences: int = None,
            allele_sequence_id: int = None,
            nb_of_samples: int = None,
    ):
        self.chrom = chrom
        self.sample = sample
        self.pos = pos
        self.ref_length = ref_length
        self.interval = interval
        self.svtype = svtype
        self.gt_conf = gt_conf
        self.coverage = coverage
        self.pangenome_variation_id = pangenome_variation_id
        self.number_of_alleles = number_of_alleles
        self.allele_id = allele_id
        self.number_of_different_allele_sequences = number_of_different_allele_sequences
        self.allele_sequence_id = allele_sequence_id
        self.nb_of_samples = nb_of_samples

    def __eq__(self, other: "ProbeHeader") -> bool:
        return (
                self.chrom == other.chrom
                and self.sample == other.sample
                and self.pos == other.pos
                and self.ref_length == other.ref_length
                and self.interval == other.interval
                and self.svtype == other.svtype
                and self.gt_conf == other.gt_conf
                and self.coverage == other.coverage
                and self.pangenome_variation_id == other.pangenome_variation_id
                and self.allele_id == other.allele_id
                and self.allele_sequence_id == other.allele_sequence_id
                and self.nb_of_samples == other.nb_of_samples
        )

    def __str__(self) -> str:
        list_of_key_values_to_add = [f"{k.upper()}={str(v)}" for k, v in vars(self).items() if v is not None]
        contents = DELIM.join(list_of_key_values_to_add)
        if not contents:
            return ""
        return f">{contents}{DELIM}"

    @staticmethod
    def from_string(string: str) -> "ProbeHeader":
        def parse_field_from_header(
                field: str, header: str, return_type: Type, value_to_return_if_not_found, delim: str = DELIM,
        ):
            regex = re.compile(f"{field}=(.+?){delim}")
            match = regex.search(header)
            if match:
                return return_type(match.group(1))
            else:
                return value_to_return_if_not_found

        chrom = parse_field_from_header("CHROM", string, str, None)
        sample = parse_field_from_header("SAMPLE", string, str, None)
        pos = parse_field_from_header("POS", string, int, None)
        ref_length = parse_field_from_header("REF_LENGTH", string, int, None)
        svtype = parse_field_from_header("SVTYPE", string, str, None)
        gt_conf = parse_field_from_header("GT_CONF", string, float, None)
        interval = ProbeInterval.from_string(
            parse_field_from_header("INTERVAL", string, str, None)
        )
        coverage = parse_field_from_header("COVERAGE", string, float, None)
        pangenome_variation_id = parse_field_from_header("PANGENOME_VARIATION_ID", string, int, None)
        number_of_alleles = parse_field_from_header("NUMBER_OF_ALLELES", string, int, None)
        allele_id = parse_field_from_header("ALLELE_ID", string, int, None)
        number_of_different_allele_sequences = parse_field_from_header("NUMBER_OF_DIFFERENT_ALLELE_SEQUENCES", string,
                                                                       int, None)
        allele_sequence_id = parse_field_from_header("ALLELE_SEQUENCE_ID", string, int, None)
        nb_of_samples = parse_field_from_header("NB_OF_SAMPLES", string, int, None)

        return ProbeHeader(
            sample=sample,
            chrom=chrom,
            pos=pos,
            ref_length=ref_length,
            interval=interval,
            svtype=svtype,
            gt_conf=gt_conf,
            coverage=coverage,
            pangenome_variation_id=pangenome_variation_id,
            number_of_alleles=number_of_alleles,
            allele_id=allele_id,
            number_of_different_allele_sequences=number_of_different_allele_sequences,
            allele_sequence_id=allele_sequence_id,
            nb_of_samples=nb_of_samples
        )


class Probe:
    def __init__(self, header: ProbeHeader = ProbeHeader(), full_sequence: str = ""):
        self.header = header
        self.full_sequence = full_sequence

    def __eq__(self, other: "Probe"):
        return self.header == other.header and self.full_sequence == other.full_sequence

    def __str__(self) -> str:
        header = str(self.header)

        if not header:
            return ""

        return f"{header}\n{self.full_sequence}"

    @property
    def left_flank(self) -> str:
        end = self.header.interval.start

        return self.full_sequence[:end]

    @property
    def right_flank(self) -> str:
        start = self.header.interval.end

        return self.full_sequence[start:]

    @property
    def core_sequence(self) -> str:
        return self.full_sequence[slice(*self.interval)]

    @property
    def interval(self) -> ProbeInterval:
        return self.header.interval

    def get_interval_or_default_interval_if_none(self) -> ProbeInterval:
        return self.interval if self.interval is not None else ProbeInterval()

    @property
    def gt_conf(self) -> float:
        return self.header.gt_conf

    @property
    def chrom(self) -> str:
        return self.header.chrom

    @property
    def pos(self) -> int:
        return self.header.pos

    @property
    def is_deletion(self) -> bool:
        return len(self.interval) == 0

    @staticmethod
    def from_string(string: str) -> "Probe":
        if not string:
            return Probe()

        fields = string.split("\n")
        header = ProbeHeader.from_string(fields[0])
        probe = Probe(header=header)

        if len(fields) > 1 and fields[1]:
            probe.full_sequence = fields[1].rstrip()

        return probe
