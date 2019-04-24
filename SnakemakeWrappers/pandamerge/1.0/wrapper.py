"""Snakemake wrapper for Pandas Merge"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2018, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import pandas as pd
import os.path as op

from snakemake.utils import makedirs

makedirs(op.dirname(snakemake.log[0]))

# Activate stdout logging
with open(snakemake.log[0], "w") as log:
    # Each key in the output dict is a column name
    for count in snakemake.output.keys():
        merged = None
        # Each key in the input dict is a sample name
        for sample, path in snakemake.input.items():
            log.write("\nWorking on {} -> {}".format(sample, path))

            # Load dataset
            data = pd.read_csv(
                path,
                sep="\t",
                index_col=0,
                header=0
            )

            # Reduce dataset
            data = data[[count]]
            data.columns = [sample]

            # Merge datasets
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

        # Build output dir if needed
        if not op.dirname(snakemake.output[count]).exists():
            makedirs(op.dirname(snakemake.output[count]))

        # Save result
        merged.to_csv(
            snakemake.output[count],
            sep="\t",
            index=True,
            header=True
        )
