"""Snakemake wrapper for Salmon Index."""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2018, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import pandas as pd
import os.path as op

from snakemake.utils import makedirs

makedirs(op.dirname(snakemake.input[0]))

merged = None
for count in snakemake.output.keys():
    for path, sample in snakemake.input.items():
        print("Working on {} -> {}".format(sample, path))
        data = pd.read_csv(
            path,
            sep="\t",
            index_col=0,
            header=0,
            dtype={
                0: str,
                1: pd.np.float,
                2: pd.np.float,
                3: pd.np.float,
                4: pd.np.float
            },
            na_values=""
        )

        data = data[[count]]
        data.columns = sample
        try:
            merged = pd.merge(
                merged,
                data,
                left_index=True,
                right_index=True
            )
        except AttributeError:
            merged_frame = data
        del data

    merged.to_csv(
        snakemake.output[count],
        sep="\t",
        index=True,
        header=True
    )
