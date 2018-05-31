#!/usr/bin/python3
# conding: utf-8

"""Plot a violin plot of the variation of length of selected genes"""

import begin
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import traceback


@begin.convert(paths=str, gencode=begin.utils.tobool)
def readSalmon(*paths: "Multiple paths to Salmon files",
               gencode: "Remove dots in target id"=False):
    """Return a DataFrame from salmon files"""
    frames = None
    for path in paths:
        logging.debug("Reading Salmon file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        data = data[["EffectiveLength"]]
        data.columns = [path]
        try:
            frames = pd.merge(frames, data, left_index=True, right_index=True)
        except ValueError:
            frames = data

    frames = pd.DataFrame(frames.stack())
    frames.reset_index(inplace=True)
    frames.columns = ["target_id", "Sample", "Length"]

    if gencode:
        frames["target_id"] = frames["target_id"].str.extract(
            '(.+).',
            expand=True
        )
        frames["target_id"] = frames["target_id"].str.rstrip('.')
    frames.set_index("target_id", inplace=True)

    return frames


@begin.convert(path=str, pvalue=float, qvalue=float,
               gencode=begin.utils.tobool)
def readSleuth(path: "Path to a Sleuth file",
               pvalue: "P-value cutoff"=0.01,
               qvalue: "Q-Value cutoff"=0.05,
               gencode: "Remove dots in target id"=False):
    """Return a DataFrame from a Sleuth file"""
    logging.debug("Reading Sleuth file: %s" % path)

    # Load DataFrame
    data = pd.read_csv(path, sep="\t", index_col=0, header=0)

    # Reduce dataframe
    data = data[data["pval"] <= pvalue]
    data = data[data["qval"] <= qvalue]

    data = pd.DataFrame(data[["b"]])
    data.columns = ["Log2FC"]
    data["Regulation"] = np.where(data["Log2FC"] >= 0,
                                  "Up-Regulated",
                                  "Down-Regulated")

    if gencode:
        data.reset_index(inplace=True)
        data.columns = ["target_id", "Log2FC", "Regulation"]
        data["target_id"] = data["target_id"].str.extract('(.+).', expand=True)
        data["target_id"] = data["target_id"].str.rstrip('.')
        data.set_index("target_id", inplace=True)

    return data


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(salmon_paths=str, sleuth_path=str, title=str, output=str,
               pvalue=float, qvalue=float, gencode=begin.utils.tobool,
               show=begin.utils.tobool)
def main(sleuth_path: "Path to a Sleuth file",
         *salmon_paths: "Path to multiple salmon files",
         title: "Title of the graph"="Variation of the length of "
                                     "the given genes",
         output: "Path to the output file"="Length_violin.png",
         pvalue: "P-value cutoff"=0.01,
         qvalue: "Q-Value cutoff"=0.05,
         gencode: "Remove dots in transcripts names"=False,
         show: "Should the graph be displayed instead of saved ?"=False):
    """Plot a violin plot of the variation of length of selected genes"""
    # Load datasetss
    frame = readSalmon(*salmon_paths, gencode=gencode)
    logging.debug(frame.head())
    sleuth = readSleuth(sleuth_path, pvalue, qvalue, gencode)
    logging.debug(sleuth.head())
    data = pd.merge(frame, sleuth, left_index=True, right_index=True)
    logging.debug(data.head())

    # Prepare plot
    logging.debug("Plotting to %s" % ("screen" if show else ": %s" % output))
    sns.set(style="whitegrid", palette="pastel", color_codes=True)
    colors = {"Up-Regulated": "b", "Down-Regulated": "y"}

    # Draw plot
    ax = sns.violinplot(x="Sample", y="Length", hue="Regulation", data=data,
                        split=True, palette=colors, inner="quart")
    ax.set_title("%s\n(p-value = %s, q-value = %s)" % (
        title,
        str(pvalue),
        str(qvalue)
    ))
    sns.despine(left=True)
    plt.xticks(rotation=45)
    plt.legend(frameon=True).get_frame().set_facecolor('white')

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
