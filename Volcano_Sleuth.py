#!/usr/bin/python3
# conding: utf-8

"""Plot a volcanoplot from a sleuth file

log2 FC
-log2 p-val

Unlincence terms of use:
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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import traceback


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(sleuth_path=str, title=str, output=str,
               qval=begin.utils.tobool, show=begin.utils.tobool)
def main(sleuth_path: "Path to a sleuth file",
         title: "Title of the graph"="Volcano plot",
         output: "Path to the output file"="Volcanoplot.png",
         qval: "Use Q-Value instead of P-Value in the volcano plot"=False,
         plot: "Should the graph be displayed instead of saved ?"=False,
         show: "Should the graph be displayed instead of saved ?"=False):
    """Plot a volcanoplot from a sleuth file"""
    logging.debug("Reading Sleuth file: %s" % sleuth_path)

    # Load DataFrame
    data = pd.read_csv(sleuth_path, sep="\t", index_col=0, header=0)

    # Reduce dataset
    col = "qval" if qval else "pval"
    col_name = "-Log10(%s)" % ("Q-Value" if qval else "P-Value")
    data = data[[col, "b"]]
    data[col] = -np.log10(data[col])
    data.columns = [col_name, "Log2(FC)"]
    logging.debug(data.head())

    logging.debug("Building graph")

    # Style
    sns.set(style="ticks", palette="pastel", color_codes=True)

    # Content
    g = sns.JointGrid(x="Log2(FC)", y=col_name, data=data)
    g = g.plot_joint(plt.scatter, color="black", edgecolor="black")
    try:
        _ = g.ax_marg_x.hist(data["Log2(FC)"], color="b", bins=30)
        _ = g.ax_marg_y.hist(data[col_name], color="r",
                             orientation="horizontal", bins=30)
    except ValueError:
        _ = g.ax_marg_x.hist(data["Log2(FC)"].dropna(), color="b", bins=30)
        _ = g.ax_marg_y.hist(data[col_name].dropna(), color="r", bins=30,
                             orientation="horizontal")

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
