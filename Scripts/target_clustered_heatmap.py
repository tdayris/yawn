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
               output=str, max_target=int, condition_array=str,
               show=begin.utils.tobool, yellow_blue=begin.utils.tobool)
def main(sleuth: "Path to a sleuth tsv file",
         expression: "Path to an expression file",
         *condition_array: "Space separated list of condition per sample",
         qval: "Q-value threashold" = 0.01,
         max_target: "Maximum number of targets" = 100,
         output: "Path to output file" = "Target_clustermap.png",
         blue: "Instead of blue->red colors" = False,
         show: "Should the graph be displayed instead of saved ?" = False):
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

    expression = expression[
        list(expression.select_dtypes(include=[np.number]).columns.values)
    ]

    expression = pd.merge(
        sleuth[["ext_gene"]],
        expression,
        left_index=True,
        right_index=True
    )

    expression.set_index("ext_gene", inplace=True)
    logging.debug(expression.head())

    # Horrible bunch of code to perform multilevel indexing
    try:
        expression = expression.T
        expression["Conditions"] = condition_array
    except ValueError:
        print(expression.index.tolist(), condition_array)
        raise

    expression = (expression.reset_index()
                            .set_index(["Conditions", "index"])
                            .T)

    # Create a categorical palette to identify the networks
    condition_lut = dict(zip(
        map(str, set(condition_array)),  # Unique conditions
        sns.husl_palette(len(set(condition_array)), s=.45)  # Color palette
    ))

    # Convert the palette to vectors
    condition_colors = (pd.Series(condition_array, index=expression.columns)
                          .map(condition_lut))

    yl = sns.diverging_palette(
        h_neg=70,
        h_pos=240,
        sep=40,
        as_cmap=True
    )

    # Draw the full plot
    ax = sns.clustermap(
        expression,
        z_score=0,
        # robust=True,
        xticklabels=True,
        yticklabels=True,
        cmap=("vlag" if not blue else "blues"),
        linewidths=.75,
        col_colors=condition_colors,
        figsize=(
            (15 if len(sleuth.columns.tolist()) <= 50 else 30),
            (15 if len(sleuth.index.tolist()) < 50 else 30)
        )
    )

    plt.setp(
        ax.ax_heatmap.xaxis.get_majorticklabels(),
        rotation=90
    )

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches="tight")
