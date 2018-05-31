#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Intention:
- Plot overlapping densities

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

import argparse
import begin
import logging
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import traceback


def label(x, color, label):
    """
    A simple function to label the plot in axes coordinates
    """
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


@begin.start(formatter_class=argparse.RawTextHelpFormatter)
@begin.logging
@begin.tracebacks
@begin.convert(distributions=str, output=str, show=begin.utils.tobool)
def main(distributions: "Path to the tsv file containing distributions",
         output: "Path to output file"="Overlapping_densities.png",
         show: "Should the graph be displayed of saved directly ?"=False):
    """
    Plot overlapping densities -- aka. Joy plot !
    """
    # Loading datasets
    data = (pd.read_csv(distributions, sep="\t", header=0, index_col=0)
              .stack()
              .reset_index())
    data.columns = ["target_id", "Samples", "Counts"]
    logging.debug("Data loaded")
    logging.debug(data.head())

    # Set environment
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)

    # Build graph
    g = sns.FacetGrid(
        data, row="Samples", hue="Samples",  # To change colors AND plot values
        aspect=15, size=.5, palette=pal
    )

    # Deal with densities and their shadows
    g.map(
        sns.kdeplot, "Counts", clip_on=False,
        shade=True, alpha=1, lw=1.5, bw=.2
    )
    g.map(
        sns.kdeplot, "Counts", clip_on=False,
        color="w", lw=2, bw=.2
    )
    g.map(plt.axhline, y=0, lw=2, clip_on=False)
    logging.debug("Data mapped")

    # Label axes
    g.map(label, "Counts")
    g.fig.subplots_adjust(hspace=-.25)
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    if show:
        logging.debug("Plotting graph ...")
        plt.show()
    else:
        logging.debug("Saving graph ...")
        plt.savefig(output, bbox_inches="tight")
