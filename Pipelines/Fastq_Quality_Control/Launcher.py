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

import argparse
import dataclasses
import os
import snakemake
import typing
import yaml


@dataclasses.dataclass(init=True, eq=False, order=False,
                       unsafe_hash=False, frozen=False)
class Fastq_Quality_Control(object):
    """
    This pipeline handles the whole quality control process

    It includes:
    - FastQC
    - FastqScreen
    - MultiQC
    """
    # Dedicated arguments
    fastq: typing.List[str] = dataclasses.field(
        init=True,
        repr=True,
        metadata={
            "help": "Space separated list of fastq files",
            "metavar": "PATH",
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
        default="/home/tdayris/Documents/Developpement/yawn/SnakemakeWrappers/cp/8.25/cold_storage.yaml",
        init=True,
        repr=True,
        metadata={
            "help": "Path to a yaml file containing the list of cold storages",
            "metavar": "PATH"
        }
    )

    # General pipeline arguments
    workdir: str = dataclasses.field(
        default="Fastq_Quality_Control_Output",
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

    # Snakemake arguments
    snakefile: str = dataclasses.field(
        init=False,
        repr=True,
        metadata={
            "help": "Path to the snakefile",
            "metavar": "PATH"
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

        # Create output directory if needed
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        self.sample = {
            os.path.basename(f): f
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
            f"cluster_config.json"
        )

        # Print self for debugging purpose
        if not self.quiet:
            print(self)

    def write(self):
        """Create the yaml config file"""
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

        # Iterate over object annotation
        for argument in cls.__annotations__.keys():

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


if __name__ == "__main__":
    args = Fastq_Quality_Control.argparser().parse_args()
    qc = Fastq_Quality_Control(**vars(args))
    qc.write()
    qc.run()
