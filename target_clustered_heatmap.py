#!/usr/bin/python3
# -*- coding: utf-8 -*-

import begin
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import traceback


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(sleuth=str, expression=str, qval=float,
               output=str, max_target=int, condition_array=str)
def main(sleuth: "Path to a sleuth tsv file",
         expression: "Path to an expression file",
         *condition_array: "Space separated list of condition per sample",
         qval: "Q-value threashold"=0.01,
         max_target: "Maximum number of targets"=100,
         output: "Path to output file"="Target_clustermap.png"):
    # Loas sleuth data
    sleuth = pd.read_csv(
        sleuth,
        sep="\t",
        header=0,
        index_col=0
    )
    sleuth = sleuth[["qval", "ext_gene"]]
    logging.debug(sleuth.head())

    sleuth = sleuth[sleuth["qval"] <= qval]
    if len(sleuth > max_target):
        sleuth = sleuth.head(n=max_target)

    # Prepare expression values
    expression = pd.read_csv(
        expression,
        sep="\t",
        header=0,
        index_col=0
    )

    expression = pd.merge(
        sleuth[["ext_gene"]],
        expression,
        left_index=True,
        right_index=True
    )

    expression.set_index("ext_gene", inplace=True)
    logging.debug(expression.head())

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=240, h_pos=10, as_cmap=True)

    # Create a categorical palette to identify the networks
    condition_lut = dict(zip(
        map(str, set(condition_array)),  # Unique conditions
        sns.husl_palette(len(set(condition_array)), s=.45)  # Color palette
    ))

    # Horrible bunch of code to perform multilevel indexing
    expression = expression.T
    expression["Conditions"] = condition_array
    expression = (expression.reset_index()
                            .set_index(["Conditions", "index"])
                            .T)

    # Convert the palette to vectors
    condition_colors = (pd.Series(condition_array, index=expression.columns)
                          .map(condition_lut))

    # Draw the full plot
    ax = sns.clustermap(
        expression,
        z_score=0,
        robust=True,
        xticklabels=True,
        yticklabels=True,
        cmap="vlag",
        linewidths=.75,
        col_colors=condition_colors,
        figsize=(min(len(sleuth), max_target), 15)
    )

    plt.setp(
        ax.ax_heatmap.xaxis.get_majorticklabels(),
        rotation=90
    )

    plt.show()
