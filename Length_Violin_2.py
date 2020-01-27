#!/usr/bin/python3
# conding: utf-8

"""
Plot a violin plot of the variation of length of differentially expressed
transcripts across multiple differential transcript analysis
"""

import begin
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import traceback

Frame = pd.DataFrame


@begin.convert(path=str, pvalue=float)
def readSleuths2(*paths: "Path to multiple Sleuth files",
                 qvalue: "Q-Value cutoff"=0.05) -> Frame:
    """Return a DataFrame from multiple Sleuth files"""
    frames = None
    for path in paths:
        logging.debug("Reading Sleuth file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        data["Regulation"] = [
            ("Up-Regulated" if x > 0 else "Down-Regulated")
            if y <= qvalue else "Equally-Expressed"
            for x, y in zip(data["b"], data["qval"])
        ]
        data = data[["Regulation"]]
        data.columns = [path]
        try:
            frames = pd.merge(
                frames,
                data,
                left_index=True,
                right_index=True,
                how='outer'
            )
        except ValueError:
            frames = data

    frames = pd.DataFrame(frames.stack()).reset_index()
    frames.columns = ["target_id", "Design", "Regulation"]
    frames.set_index("target_id", inplace=True)

    return frames


@begin.start
@begin.logging
@begin.tracebacks
def main(length_path: "Path to transcript length file",
         *sleuth_path: "Path to multiple sleuth files",
         title: "Title of the graph"="Variation of the size of transcripts"
                                     " across experimental designs",
         output: "Path to the output file"="Length_across_Design.png",
         qvalue: "Q-Value cutoff"=0.05,
         show: "Plot the graph on screen instead of saving it"=True) -> None:
    """
    Plot a violin plot of the variation of length of differentially expressed
    transcripts across multiple differential transcript analysis
    """
    # Geather data
    sleuth = readSleuths2(*sleuth_path)
    lengths = pd.read_csv(length_path, sep="\t", header=0, index_col=0)
    #
    logging.debug(sleuth.head())
    logging.debug(lengths.head())
    #
    data = pd.merge(sleuth, lengths, left_index=True, right_index=True)
    logging.debug(data.head())

    g = sns.boxplot(
        x="Design",
        y="Length",
        hue="Regulation",
        data=data,
        palette="muted"
    )

    plt.xticks(rotation=10)
    plt.yticks(np.arange(0, data["Length"].max(), 1000))
    plt.title(title)
    sns.despine(left=True, trim=True)

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
