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

__all__ = ["LauncherSnakemake"]
__version__ = "0.1.1"

import argparse
import dataclasses
import itertools
import os
import snakemake
import typing
import yaml


@dataclasses.dataclass(init=True, eq=False, order=False,
                       unsafe_hash=False, frozen=False)
class LauncherSnakemake(object):
    """
    This object is a metaclass which will help further pipelines to be made
    Just inherit from this object and the following functions will be bought
    to your pipeline:

    - Automatic argument parser from class annotation
    - Multiple Snakemake parameters + call of snakemake itself
    - Default singularity image
    """
    # General snakemake parameters without default argument
    snakefile: str = dataclasses.field(
        init=False,
        repr=True,
        metadata={
            "help": "Path to the snakefile",
            "metavar": "PATH"
        }
    )

    # General pipeline arguments with default arguments
    workdir: str = dataclasses.field(
        default="Snakemake_output",
        init=True,
        repr=True,
        metadata={
            "help": "Path to the output directory",
            "metavar": "PATH",
            "short": "w"
        }
    )
    threads: int = dataclasses.field(
        default=1,
        init=True,
        repr=True,
        metadata={
            "help": "Maximum number of threads used by a process",
            "short": "t"
        }
    )
    quiet: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Reduce verbosity",
            "action": "store_true"
        }
    )
    jobs: int = dataclasses.field(
        default=20,
        init=True,
        repr=True,
        metadata={
            "help": "Maximum number of jobs submitted"
                    " (ignored in local execution)",
            "short": "j"
        }
    )
    dryrun: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Perform a dry-run, without any execution",
            "action": "store_true",
            "short": "n"
        }
    )
    force: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Force all rules to be executed. This overwrites possible"
                    "existing results.",
            "action": "store_true",
            "short": "f"
        }
    )
    keepgoing: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Continue to go through the execution tree, eventhough "
                    "some jobs have failed",
            "action": "store_true",
            "short": "k"
        }
    )
    dag: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Pring DAG of the execution tree",
            "action": "store_true"
        }
    )
    graph: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Pring the gule graph corresponding to the execution tree",
            "action": "store_true"
        }
    )
    keeptmp: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Do not delete temporary files",
            "action": "store_true"
        }
    )
    singularity: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Use singularity and run the whole pipeline in containers",
            "action": "store_true"
        }
    )
    noconda: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Do not use conda to execute the pipeline",
            "action": "store_true"
        }
    )
    restart: int = dataclasses.field(
        default=3,
        init=True,
        repr=True,
        metadata={
            "help": "Number of restart time for failed jobs"
        }
    )
    summary: bool = dataclasses.field(
        default=False,
        init=True,
        repr=True,
        metadata={
            "help": "Print a detailed summary of the whole execution",
            "action": "store_true"
        }
    )
    singularity_docker_image: str = dataclasses.field(
        default="docker://continuumio/miniconda3:4.4.10",
        repr=True,
        init=True,
        compare=False,
        metadata={
            "help": "Address of the docker/singularity image (it MUST "
                    "contain conda or else, you have to use --noconda)"
        }
    )

    def __post_init__(self):
        """
        Post init instructions

        Remember to set a snakefile: self.snakefile = ...
        """
        self.snakefile = ""

    def write(self) -> None:
        """Create the yaml config file"""
        # Create output directory if needed
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # No control over config file overwriting
        with open(self.config, "w") as config_file:
            data = yaml.dump(
                dataclasses.asdict(self),
                default_flow_style=False
            )
            config_file.write(data)

    def run(self):
        """Launch the snakemake through its API"""
        # Snakemake already moves into output directory
        # I will try to remove it
        os.chdir(self.workdir)

        snakemake.snakemake(
            snakefile=self.snakefile,
            cores=self.threads,
            nodes=self.jobs,
            local_cores=min(2, self.threads),
            configfile=self.config,
            workdir=self.workdir,
            dryrun=self.dryrun,
            printshellcmds=True,
            cluster_config=self.cluster_config,
            printreason=True,
            keepgoing=self.keepgoing,
            jobname="{rulename}.snakejob.{jobid}",
            latency_wait=60,
            snakemakepath="/cm/shared/bioinfo/python/3.5.2/bin/snakemake",
            printdag=self.dag,
            quiet=self.quiet,
            notemp=self.keeptmp,
            forceall=self.force,
            use_singularity=self.singularity,
            printrulegraph=self.graph,
            use_conda=(not self.noconda),
            summary=self.summary,
            restart_times=self.restart
        )

    @classmethod
    def argparser(cls) -> argparse.ArgumentParser:
        """
        Return the command line parser for the pipeline
        """
        # Build the main parser
        parser_pipeline = argparse.ArgumentParser(
            description=cls.__doc__,
            epilog="All tools used within this pipeline belong to their "
                   "respective authors.",
            formatter_class=argparse.RawTextHelpFormatter
        )

        args_iter = itertools.chain(
            cls.__annotations__.keys(),
            cls.__base__.__annotations__.keys()
        )
        # Iterate over object annotation
        for argument in args_iter:

            # If an argument starts with "_", then
            # it is "private". We shall ignore it.
            public = not argument.startswith("_")

            # Get the corresponding dataclass field
            value = cls.__dataclass_fields__[argument]

            # If an argument is not in the __init__ function
            # Then it does not belong to be in the command
            # line parser.
            inited = value.init is True

            if public and inited:
                # dataclasses._MISSING_TYPE is the default
                # value, in a dataclasses field.
                required = isinstance(
                    value.default,
                    dataclasses._MISSING_TYPE
                )
                if argument == "fastq":
                    print(value)

                # The help message written in the dataclasses field
                help = f"{value.metadata['help']}"

                # Optional argument starts with two dashes (--)
                args = [f"{argument}"]
                if not required:
                    args = [f"--{argument}"]
                    help += f" (default: {value.default})"
                    if "short" in value.metadata.keys():
                        args.append(f"-{value.metadata['short']}")

                # The argument parser will require multiple
                # arguments, we will store them in this variables
                kwargs = {
                    "help": help,
                    "default": value.default
                }

                # Now, according to the type of the dataclasses field,
                # the argument parser will be filled differently
                basic_types = value.type in [str, int, float]
                is_bool = value.type == bool
                if is_bool:
                    # Type must not be set in booleans
                    # Default value must not be set in argparse boolean
                    kwargs["action"] = value.metadata["action"]

                if basic_types:
                    # Only for basic types
                    kwargs["type"] = value.type

                if "metavar" in value.metadata.keys():
                    kwargs["metavar"] = value.metadata["metavar"]

                if "choices" in value.metadata.keys():
                    kwargs["choices"] = value.metadata["choices"]

                kwargs["help"] = help

                if "nargs" in value.metadata.keys():
                    kwargs["nargs"] = value.metadata["nargs"]

                # Eventually, the argument must be added
                parser_pipeline.add_argument(*args, **kwargs)

        return parser_pipeline
