"""Snakemake wrapper for Salmon Index."""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2018, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import pandas as pd
import os.path as op

from snakemake.utils import makedirs

makedirs(op.dirname(snakemake.input[0]))
makedirs(op.dirname(snakemake.log[0]))

with open(snakemake.log[0], "w") as log:
    for count in snakemake.output.keys():
        merged = None
        for sample, path in snakemake.input.items():
            log.write("\nWorking on {} -> {}".format(sample, path))
            data = pd.read_csv(
                path,
                sep="\t",
                index_col=0,
                header=0
            )

            data = data[[count]]
            data.columns = [sample]
            try:
                merged = pd.merge(
                    merged,
                    data,
                    left_index=True,
                    right_index=True
                )
            except TypeError:
                merged = data
            log.write(str(merged.head()))
            log.write("\n")
            del data

        merged.to_csv(
            snakemake.output[count],
            sep="\t",
            index=True,
            header=True
        )
