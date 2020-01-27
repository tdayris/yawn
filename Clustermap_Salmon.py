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
import os.path as op
import pandas as pd
import seaborn as sns
import traceback

# from matplotlib import interactive
# interactive(True)


@begin.convert(paths=str, prefix=str, use_dir_name=begin.utils.tobool)
def readSalmon(*paths: "Multiple paths to Salmon files",
               prefix: "Remove prefix from sample names"="SAQuant_",
               use_dir_name: "Use dirname instead of file name "
                             "for samples"=False):
    """Return a DataFrame from salmon files"""
    frames = None
    for path in paths:
        logging.debug("Reading Salmon file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        name = path
        if use_dir_name:
            name = op.basename(op.dirname(path))
            name = name.lstrip(prefix)
        data = data[["TPM"]]
        data.columns = [name]
        try:
            frames = pd.merge(frames, data, left_index=True, right_index=True)
        except ValueError:
            frames = data

    return frames


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(salmon_paths=str, output=str, prefix=str, label_rotation=int,
               show=begin.utils.tobool, use_dir_name=begin.utils.tobool)
def main(*salmon_paths: "Path to multiple Salmon files",
         output: "Path to the output file"="Clustered_Heatmap.png",
         show: "Should the graph be displayed instead of saved ?"=False,
         prefix: "Remove prefix from sample names"="SAQuant_",
         label_rotation: "The label rotation used to avoid superposition"=20,
         count_table: "The input path is a unique count table,"
                      " not a serie of salmon files"=False,
         use_dir_name: "Use dirname instead of file name for samples"=False):
    """Plot a clustered heatmap of multiple sample"""

    # Load datasets
    frame = (
        pd.read_csv(*salmon_paths, sep="\t", header=0, index_col=0)
        if count_table else
        readSalmon(*salmon_paths, use_dir_name=use_dir_name, prefix=prefix)
    )

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=240, h_pos=10, as_cmap=True)

    # Draw the full plot
    ax = sns.clustermap(frame.corr(), linewidths=0.5,
                        figsize=(13, 13), cmap=cmap)

    plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(ax.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches="tight")
