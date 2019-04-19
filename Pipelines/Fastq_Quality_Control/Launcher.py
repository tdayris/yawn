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

__all__ = ["Fastq_Quality_Control"]
__version__ = "0.1.1"

# General importations
import argparse
import dataclasses
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
class Fastq_Quality_Control(LauncherSnakemake):
    """
    This pipeline handles the whole quality control process

    It includes:
    - FastQC
    - FastqScreen
    - MultiQC
    """
    # Dedicated arguments
    fastq: typing.List[str] = dataclasses.field(
        default=dataclasses._MISSING_TYPE(),
        init=True,
        repr=True,
        metadata={
            "help": "Space separated list of fastq files",
            "metavar": "FASTQ",
            "nargs": "+"
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
        if self.fastq[0].endswith(("tar.gz", "tar.bz2")):
            self.fq_ext = ".".join(self.fastq[0].split(".")[-3:])
        elif self.fastq[0].endswith(("gz", "bz2")):
            self.fq_ext = ".".join(self.fastq[0].split(".")[-2:])
        else:
            self.fq_ext = self.fastq[0].split(".")[-1]

        if not os.path.isabs(self.workdir):
            self.workdir = os.path.realpath(self.workdir)

        self.sample = {
            os.path.basename(f): os.path.realpath(f)
            for f in self.fastq
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
    print(Fastq_Quality_Control.__annotations__)
    args = Fastq_Quality_Control.argparser().parse_args()
    qc = Fastq_Quality_Control(**vars(args))
    qc.write()
    qc.run()
