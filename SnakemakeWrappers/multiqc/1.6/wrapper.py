#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

__author__ = "Thibault Dayris"
__copyright__ = "Fish_n_CHIP 2019"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "GPLV3"

import os

from snakemake.shell import shell

# Format input files
multiqc_input = ""
if isinstance(snakemake.input, list):
    multiqc_input = " ".join(snakemake.input)
elif isinstance(snakemake.input, dict):
    multiqc_input = " ".join(snakemake.input.values())
else:
    multiqc_input = snakemake.input

# Get output directory
output_directory = os.path.dirname(snakemake.output[0])
output_name = os.path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell_cmd = (
    "multiqc "                          # Tool
    "{snakemake.params.extra} "         # Extra parameters
    "--force "                          # Overwrite existing reports
    "--outdir {output_directory} "      # Output directory
    "--filename {output_name} "         # Output file name
    "{multiqc_input} "                  # Space separated list of files
    "{log}"                             # Logging
)

shell(shell_cmd)
