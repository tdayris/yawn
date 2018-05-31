#!/usr/bin/python3
# conding: utf-8

"""Plot a pairwise scatterplot of multiple samples"""

import begin
import logging
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import traceback


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(salmon_paths=str, output=str, show=begin.utils.tobool)
def main(*salmon_paths: "Path to salmon files",
         output: "Path to the output file"="Pairwise_scatterplot.png",
         show: "Should the graph be displayed instead of saved ?"=False):
    """Plot a pairwise scatterplot of multiple samples"""
    # Load datasets
    frames = None
    for path in salmon_paths:
        logging.debug("Reading Salmon file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        data = data[["TPM"]]
        data.columns = [path]
        try:
            frames = pd.merge(frames, data, left_index=True, right_index=True)
        except ValueError:
            frames = data

    # Build graph
    logging.debug("Building graph")
    sns.set(style="ticks", color_codes=True)
    g = sns.pairplot(frames,
                     diag_kind="kde",
                     diag_kws=dict(shade=True))

    logging.debug("Plotting to %s" % ("standard output" if show else output))
    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
