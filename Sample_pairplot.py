#!/usr/bin/python3
# conding: utf-8

import begin
import logging
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import traceback

from pathlib import Path


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(counts=str, outplot=str, legend=str,
               dropna=begin.utils.tobool, show=begin.utils.tobool)
def main(counts: "Path to a counts tsv file",
         outplot: "Path to output graph" = "Pairplot.png"
         legend: "Column name with levels to color" = "",
         dropna: "Drop missing values from the data before plotting" = False,
         show: "Plot the graph on screen, instead of saving it" = False):
    """
    This script builds a pairplot from sample counts

    This count table shoud look like:
            Sample1     Sample2     ...     SampleN
    Gene1   0.0156      15.34       ...     47.25
    Gene2   NaN         12.58       ...     44.35
    Gene3   1.247       147.37      ...     NaN
    ...     ...         ...         ...     ...
    GeneN   1.2789      147.22      ...     45.369

    Exemple with a legend column:
            Sample1     Sample2     ...     SampleN     Factor
    Gene1   0.0156      15.34       ...     47.25       Level1
    Gene2   NaN         12.58       ...     44.35       Level1
    Gene3   1.247       147.37      ...     NaN         Level1
    ...     ...         ...         ...     ...         ...
    GeneN   1.2789      147.22      ...     45.369      LevelN
    """
    counts = Path(counts)
    if not counts.exists():
        raise FileNotFoundError("Could not find: %s" % str(counts))

    data = pd.read_csv(
        str(counts),
        sep="\t",
        header=0,
        index=0
    )

    sns.set(
        style="darkgrid"
    )

    sns.pairplot(
        data,
        hue=(None if legend == "" else legend),
        diag_kind="kde",
        palette="Blues"
    )

    if show:
        plt.tight_layout()
        plt.show()
    else:
        plt.savefig(outplot, bbox_inches='tight')
