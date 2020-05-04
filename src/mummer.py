"""This file holds wrappers for running mummer commands."""
import logging
import re
import subprocess
from pathlib import Path
from typing import List, TextIO, Tuple

import pandas as pd


class NucmerError(Exception):
    pass


class DeltaFilterError(Exception):
    pass


class ShowSnpsError(Exception):
    pass


class GetReportFromDeltaFileError(Exception):
    pass


class Nucmer:
    def __init__(
            self, reference: Path, query: Path, prefix: str = "out", extra_params: str = ""
    ):
        if not reference.is_file():
            raise NucmerError(f"Reference file {str(reference)} does not exist.")
        elif not query.is_file():
            raise NucmerError(f"Query file {str(query)} does not exist.")
        if not Path(prefix).parent.is_dir():
            raise NucmerError(f"Prefix {str(Path(prefix))} parent does not exist.")

        self.reference = str(reference)
        self.query = str(query)
        self.prefix = prefix
        self.extra_params = extra_params

    def generate_command(self) -> List[str]:
        command = ["nucmer", "--prefix", self.prefix, self.reference, self.query]

        if self.extra_params:
            command.insert(1, self.extra_params)

        return command

    def run(self) -> subprocess.CompletedProcess:
        """The output file is written to <prefix>.delta"""
        nucmer_command = self.generate_command()

        logging.info(f"Running nucmer with command:\n{' '.join(nucmer_command)}")

        return subprocess.run(
            nucmer_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )


class GetReportFromDeltaFile:
    def __init__(self, deltafile: Path):
        if not deltafile.is_file():
            raise GetReportFromDeltaFileError(f"deltafile {str(deltafile)} does not exist.")
        self.deltafile = str(deltafile)

    def generate_command(self) -> List[str]:
        command = ["dnadiff", "-d", self.deltafile, "-p", f"{self.deltafile}.dnadiff_prefix"]
        return command

    def get_report(self) -> str:
        """The output file can be found in the stdout of the returned
        CompletedProcess."""
        command = self.generate_command()
        subprocess.check_call(command)
        report_file = Path(f"{self.deltafile}.dnadiff_prefix.report")
        return report_file.read_text()

    @staticmethod
    def get_ref_and_query_aligned_bases_percentage(report) -> Tuple[float, float]:
        matches = re.search(r'AlignedBases\s+\d+\((\d*\.\d*)%\)\s+\d+\((\d*\.\d*)%\)', report)
        ref_aligned_bases_percentage = float(matches[1])
        query_aligned_bases_percentage = float(matches[2])
        return ref_aligned_bases_percentage, query_aligned_bases_percentage


class DeltaFilter:
    def __init__(self, deltafile: Path, extra_params: str = ""):
        if not deltafile.is_file():
            raise DeltaFilterError(f"deltafile {str(deltafile)} does not exist.")

        self.deltafile = str(deltafile)
        self.extra_params = extra_params

    def generate_command(self) -> List[str]:
        command = ["delta-filter", self.deltafile]

        if self.extra_params:
            command.insert(1, self.extra_params)

        return command

    def run(self) -> subprocess.CompletedProcess:
        """The output file can be found in the stdout of the returned
        CompletedProcess."""
        deltafilter_command = self.generate_command()

        logging.info(
            f"Running delta-filter with command:\n{' '.join(deltafilter_command)}"
        )

        return subprocess.run(
            deltafilter_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )


class ShowSnps:
    def __init__(
            self,
            deltafile: Path,
            context: int = 0,
            print_header: bool = True,
            indels: bool = True,
            extra_params: str = "",
    ):
        if not deltafile.is_file():
            raise ShowSnpsError(f"deltafile {str(deltafile)} does not exist.")

        self.deltafile = str(deltafile)
        self.context = context
        self.print_header = print_header
        self.indels = indels
        self.extra_params = extra_params

    def generate_command(self) -> List[str]:
        command = ["show-snps", self.deltafile]

        if not self.indels:
            command.insert(1, "-I")

        if not self.print_header:
            command.insert(1, "-H")

        if self.context > 0:
            command.insert(1, f"-x {self.context}")

        if self.extra_params:
            command.insert(1, self.extra_params)

        return command

    def run(self) -> subprocess.CompletedProcess:
        """The output file can be found in the stdout of the returned
        CompletedProcess."""
        showsnps_command = self.generate_command()

        logging.info(f"Running show-snps with command:\n{' '.join(showsnps_command)}")

        return subprocess.run(
            showsnps_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    @staticmethod
    def to_dataframe(snps: TextIO) -> "ShowSNPsDataframe":
        """Note: this method is not general. i.e it is only setup at the moment to
        parse a show-snps file where the options used were -rlTC and -x"""
        cols = {
            "ref_pos": int,  # P1,
            "ref_sub": str,  # SUB,
            "query_sub": str,  # SUB
            "query_pos": int,  # P2
            "nearest_mismatch": int,  # BUFF
            "nearest_end": int,  # DIST
            "ref_len": int,  # LEN R
            "query_len": int,  # LEN Q
            "ref_context": str,  # CTX R
            "query_context": str,  # CTX Q
            "ref_strand": int,
            "query_strand": int,
            "ref_chrom": str,
            "query_chrom": str,
        }
        names = list(cols.keys())
        return ShowSNPsDataframe(
            pd.read_csv(
                snps, sep="\t", skiprows=4, index_col=False, names=names, dtype=cols
            )
        )


class ShowSNPsDataframe(pd.DataFrame):
    @property
    def _constructor(self):
        return ShowSNPsDataframe

    def translate_to_FWD_strand(self) -> "ShowSNPsDataframe":
        def fix_position(position: int, strand_aln: int, length: int) -> int:
            if strand_aln == 1:
                return position
            else:
                return length - position + 1

        def translate_to_FWD_strand_core(line: pd.Series) -> pd.Series:
            line.ref_pos = fix_position(line.ref_pos, line.ref_strand, line.ref_len)
            line.ref_strand = 1
            line.query_pos = fix_position(
                line.query_pos, line.query_strand, line.query_len
            )
            line.query_strand = 1
            return line

        return self.apply(translate_to_FWD_strand_core, axis=1)

    def make_pos_zero_based(self) -> "ShowSNPsDataframe":
        df = self.copy(deep=True)
        if df.empty:
            return df
        df["ref_pos"] = df["ref_pos"] - 1
        df["query_pos"] = df["query_pos"] - 1
        return df
