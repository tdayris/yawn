#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for FastQC_0_11_7"""

__author__ = ["Thibault Dayris"]
__copyright__ = "Fish_N_CHIP 2017"
__email__ = ["thibault.dayris@gustaveroussy.fr"]
__license__ = "GPLV3"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell_cmd = (
    "fastqc "                                 # Tool
    "{snakemake.params.extra} "               # Extra parameters for FastQC
    "--outdir {snakemake.params.tempdir} "    # Output directory
    "--threads {snakemake.threads} "          # Number of threads used
    "{snakemake.input[0]} "                   # Path to fastq file
    "{log}"                                   # Logging
)

shell(shell_cmd)
