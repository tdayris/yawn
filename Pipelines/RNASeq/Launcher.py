#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

__author__ = "Thibault Dayris"
__copyright__ = "Fish_n_CHIP 2019"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "GPLV3"

__all__ = ["RNASeq"]
__version__ = "0.1.0"


# General importations
import argparse
import dataclasses
# import itertools
import os
import typing
import sys

try:
    sys.path.append(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    )
    from LauncherSnakemake import LauncherSnakemake
except ModuleNotFoundError:
    raise


@dataclasses.dataclass(init=True, eq=False, order=False,
                       unsafe_hash=False, frozen=False)
class RNASeq(LauncherSnakemake):
    """
    This pipeline handles the whole RNA-Seq analysis

    It includes:
    - Salmon
    """
    fastq_r1: typing.List[str] = dataclasses.field(
        default="",
        init=True,
        repr=False,
        metadata={
            "help": "Space separated list of fastq files",
            "metavar": "FASTQ",
            "nargs": "+",
            "short": "r1"
        }
    )
    fastq_r2: typing.List[str] = dataclasses.field(
        default="",
        init=True,
        repr=False,
        metadata={
            "help": "Space separated list of fastq files",
            "metavar": "FASTQ",
            "nargs": "+",
            "short": "r2"
        }
    )
    conditions: typing.List[str] = dataclasses.field(
        default="",
        init=True,
        repr=False,
        metadata={
            "help": "Space separated list of conditions,"
                    " corresponding to the fastq files"
        }
    )
    fasta_transcriptome: str = dataclasses.field(
        default=dataclasses._MISSING_TYPE(),
        init=True,
        repr=False,
        metadata={
            "help": "Path to a fasta-formatted transcriptome sequence",
            "metavar": "TRANSCRIPTOME",
        }
    )
    fasta_genome: str = dataclasses.field(
        default=dataclasses._MISSING_TYPE(),
        init=True,
        repr=False,
        metadata={
            "help": "Path to a fasta-formatted genome sequence",
            "metavar": "GENOME",
        }
    )
    gtf: str = dataclasses.field(
        default=dataclasses._MISSING_TYPE(),
        init=True,
        repr=False,
        metadata={
            "help": "Path to a gtf-formatted genome annotation",
            "metavar": "GTF",
        }
    )
    sample: dict = dataclasses.field(
        default_factory=dict,
        init=False,
        repr=True,
        metadata={
            "help": "Link between sample names and files paths"
        }
    )
    sources: dict = dataclasses.field(
        default_factory=dict,
        init=False,
        repr=True,
        metadata={
            "help": "Link between genome information, and files paths"
        }
    )
    pair: dict = dataclasses.field(
        default_factory=dict,
        init=False,
        repr=True,
        metadata={
            "help": "Link between pairs of fastq files"
        }
    )
    salmon_index_extra: str = dataclasses.field(
        default="",
        init=True,
        repr=True,
        metadata={
            "help": "Extra parameters for salmon index. "
                    "Write them between simple quotes."
        }
    )
    salmon_quant_extra: str = dataclasses.field(
        default="",
        init=True,
        repr=True,
        metadata={
            "help": "Extra parameters for salmon quant. "
                    "Write them between simple quotes."
        }
    )
    libtype: str = dataclasses.field(
        default="A",
        init=True,
        repr=True,
        metadata={
            "help": "The sequencing library, ask your genomic platform.",
            "choices": ["A", "ISF", "ISR", "IU", "MSF",
                        "MSR", "MU", "OSR", "OSF", "OU"]
        }
    )
    cold_storage: str = dataclasses.field(
        init=False,
        repr=True,
        metadata={
            "help": "Path to a yaml file containing the list of cold storages",
            "metavar": "PATH"
        }
    )

    # Hidden arguments
    fq_ext: str = dataclasses.field(
        init=False,
        repr=True,
        metadata={
            "help": "The extension of all fastq files "
                    "(used for snakemake regexp)",
        }
    )
    config: str = dataclasses.field(
        init=False,
        repr=True,
        metadata={
            "help": "Path to the yaml config file for snakemake",
            "metavar": "PATH"
        }
    )
    cluster_config: str = dataclasses.field(
        init=False,
        repr=True,
        metadata={
            "help": "Path to the yaml file containing"
                    " cluster specific information",
            "metavar": "PATH"
        }
    )

    def __post_init__(self):
        """Post init instructions"""
        # We need to recover the extention of all fastq files
        # This is required for snakemake wildcards treatment
        if self.fastq_r1[0].endswith(("tar.gz", "tar.bz2")):
            self.fq_ext = ".".join(self.fastq_r1[0].split(".")[-3:])
        elif self.fastq_r1[0].endswith(("gz", "bz2")):
            self.fq_ext = ".".join(self.fastq_r1[0].split(".")[-2:])
        else:
            self.fq_ext = self.fastq_r1[0].split(".")[-1]

        if not os.path.isabs(self.workdir):
            self.workdir = os.path.realpath(self.workdir)

        self.sample = {
            os.path.basename(f): os.path.realpath(f)
            for f in self.fastq_r1 + self.fastq_r2
        }

        self.sources = {
            os.path.basename(f): os.path.realpath(f)
            for f in [self.gtf, self.fasta_transcriptome, self.fasta_genome]
        }

        self.confitions = {
            sample: condition
            for sample, condition in zip(self.sample.keys(), self.conditions)
        }

        for up, down in zip(self.fastq_r1, self.fastq_r2):
            sample_id = os.path.commonprefix([
                os.path.basename(up),
                os.path.basename(down)
            ])
            self.pair[sample_id] = {
                "r1": os.path.basename(up),
                "r2": os.path.basename(down)
            }

        # Gather required paths
        self.config = os.path.join(
            self.workdir,
            f"snakemake_config_{self.__class__.__name__}.yaml"
        )
        self.snakefile = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "Snakefile"
        )
        self.cluster_config = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "cluster_config.json"
        )
        self.cold_storage = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "cold_storage.yaml"
        )


if __name__ == "__main__":
    args = RNASeq.argparser().parse_args()
    rna = RNASeq(**vars(args))
    rna.write()
    rna.run()
