#!/usr/bin/python3
# conding: utf-8

"""
Plot a clustered heatmap of multiple sample

Unlincense terms of use:
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""

import begin
import logging
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import pandas as pd
import seaborn as sns
import traceback


@begin.start
@begin.logging
@begin.tracebacks
def main(normalized_table: "Path to a normalized quantification table",
         *condition_array: "Condition name per sample (space separated)",
         title: "Graph title"="Clustered Heatmap",
         output: "Path to output file"="Clustered_Heatmap.png",
         xlabel_rotation: "The X label rotation"=0,
         ylabel_rotation: "The Y label rotation"=90,
         row_cluster: "Perform clustering on rows"=True,
         col_cluster: "Perform clustering on columns"=True,
         row_condition_color: "Set colors on rows according to "
                              "condition array"=True,
         col_condition_color: "Set colors on columns according to "
                              "condition array"=True,
         show: "Should the graph be displayed instead of saved ?"=False,
         square: "Plot half of the heatmap"=False,
         robust: "use robust quantiles instead of the extreme values"=False):
    """
    Plot a clustered heatmap of multiple sample

    normalized_table should be a path to a tsv file, where first line contains
    sample names, first column contains target id.
    """
    # Dataset loading
    logging.debug("Loading %s" % normalized_table)
    data = pd.read_csv(normalized_table, sep="\t", header=0, index_col=0)
    data = data[list(data.select_dtypes(include=[np.number]).columns.values)]
    logging.debug(data.head())
    if len(data.columns.tolist()) != len(condition_array):
        raise ValueError(
            "Expected the same number of condition and samples"
            ", got different ones %i / %i" % (
                len(data.columns.tolist()),
                len(condition_array)
            )
        )

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=240, h_pos=10, as_cmap=True)

    # Create a categorical palette to identify the networks
    condition_lut = dict(zip(
        map(str, set(condition_array)),  # Unique conditions
        sns.husl_palette(len(set(condition_array)), s=.45)  # Color palette
    ))

    # Horrible bunch of code to perform multilevel indexing
    data = data.T
    data["Conditions"] = condition_array
    data = (data.reset_index()
                .set_index(["Conditions", "index"])
                .T)

    # Convert the palette to vectors
    condition_colors = (pd.Series(condition_array, index=data.columns)
                          .map(condition_lut))

    data = data.corr()
    logging.debug(data.head())

    # Draw the full plot
    ax = sns.clustermap(
        data,
        cmap=cmap,
        row_colors=(condition_colors if row_condition_color else None),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        col_colors=(condition_colors if col_condition_color else None),
        linewidths=0.5,
        figsize=(13, 13),
        square=square,
        robust=robust
    )

    plt.setp(
        ax.ax_heatmap.yaxis.get_majorticklabels(),
        rotation=xlabel_rotation)
    plt.setp(
        ax.ax_heatmap.xaxis.get_majorticklabels(),
        rotation=ylabel_rotation
    )

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches="tight")
