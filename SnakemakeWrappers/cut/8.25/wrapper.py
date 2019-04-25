#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bash cut"""

__author__ = "Thibault Dayris"
__copyright__ = "Fish_n_CHIP 2019"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "GPLV3"

import os

from snakemake.shell import shell

# Prepare logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Defining parameters
delimiter = snakemake.params.get("delimiter", "$'\t'")
columns = snakemake.params.get("columns", "3")
extra = snakemake.params.get("extra", "")

shell(
    "cut "                   # Tool
    "-f {columns} "          # Selection of columns/fields
    "-d {delimiter} "        # Selection of delimiter
    "{extra} "               # Extra parameters
    "{snakemake.input} "     # Path to input file
    "> {snakemake.output} "  # Path to output file
    "{log}"                  # Logging
)
