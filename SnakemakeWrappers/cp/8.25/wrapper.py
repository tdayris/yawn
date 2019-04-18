#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

__author__ = "Thibault Dayris"
__copyright__ = "Fish_n_CHIP 2019"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "GPLV3"

import os
import yaml

from snakemake.shell import shell

# No node can access a cold storage
# these files must be copied. However,
# any file else where should be symlinked!
with open(snakemake.params.cold_storage, "r") as cold_list:
    cold_storage = yaml.load(cold_list)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell_cmd = ""
for f in snakemake.input:
    symlink = any(f.startswith(cs) for cs in cold_storage)
    shell_cmd += (
        "cp "                                        # Tool
        "{snakemake.params.extra} "                  # Optionnal parameters
        f"{'' if symlink else '--symbolic-link'} "   # Symlink if needed
        f"{os.path.realpath(f)} "                    # Path to input file
        f"{os.path.realpath(snakemake.output[0])} "  # Path to output file
        "{log}; "                                    # Logging
    )

shell(shell_cmd)
