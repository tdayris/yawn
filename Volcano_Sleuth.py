#!/usr/bin/python3
# conding: utf-8

"""Plot a volcanoplot from a sleuth file

log2 FC
-log2 p-val

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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import traceback


@begin.convert(stat_threshold=float, beta_threshold=float)
def classify(betas, stats, beta_threshold, stat_threshold):
    """
    Yields a classification of events
    """
    for b, s in zip(betas, stats):
        if s <= stat_threshold:
            if abs(b) >= beta_threshold:
                yield "Differentially Expressed"
            else:
                yield "Below beta threshold"
        else:
            yield "Below significance threshold"


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(sleuth_path=str, output=str, stat_threshold=float,
               beta_threshold=float, qval=begin.utils.tobool,
               show=begin.utils.tobool)
def main(sleuth_path: "Path to a sleuth file",
         output: "Path to the output file" = "Volcanoplot.png",
         stat_threshold: "The threshold above which "
                         "targets are diff. exp." = 0.05,
         beta_threshold: "The threshold above which "
                         "fold changes are significant" = 1,
         qval: "Use P-Value instead of Q-Value in the volcano plot" = True,
         show: "Should the graph be displayed instead of saved ?" = False):
    """Plot a volcanoplot from a sleuth file"""
    logging.debug("Reading Sleuth file: %s" % sleuth_path)

    # Load DataFrame
    data = pd.read_csv(sleuth_path, sep="\t", index_col=0, header=0)

    # Reduce dataset
    col = "qval" if qval else "pval"
    col_name = "-Log10(%s)" % ("Q-Value" if qval else "P-Value")
    data = data[[col, "b"]]
    data.columns = [col_name, "Log2(beta)"]
    data["Significance"] = list(classify(
        data["Log2(beta)"], data[col_name], beta_threshold, stat_threshold
    ))
    data[col_name] = [0 if d == 0 else -np.log10(d) for d in data[col_name]]
    logging.debug(data.head())

    logging.debug("Building graph")

    # Style
    sns.set(style="ticks", palette="pastel", color_codes=True)

    # Content
    g = sns.lmplot(
        x="Log2(beta)",
        y=col_name,
        hue="Significance",
        hue_order=[
            "Below significance threshold",
            "Below beta threshold",
            "Differentially Expressed"
        ],
        data=data,
        palette="Blues",
        fit_reg=False
    )
    g.map(
        plt.axhline,
        y=-np.log10(stat_threshold),
        ls=":",
        c=".5",
        markersize=0.1
    )
    g.map(
        plt.axvline,
        x=beta_threshold,
        ls=":",
        c=".5",
        markersize=0.1
    )
    g.map(
        plt.axvline,
        x=-beta_threshold,
        ls=":",
        c=".5",
        markersize=0.1
    )

    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
