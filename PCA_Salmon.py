#!/usr/bin/python3
# conding: utf-8

"""
Plot a PCA on selected samples

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
import itertools
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.decomposition as skd
import traceback


@begin.convert(paths=str)
def readSalmon(*paths: "Multiple paths to Salmon files"):
    """Return a DataFrame from salmon files"""
    frames = None
    for path in paths:
        logging.debug("Reading Salmon file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        data = data[["TPM"]]
        data.columns = [path]
        try:
            frames = pd.merge(frames, data, left_index=True, right_index=True)
        except ValueError:
            frames = data

    return frames


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(salmon_paths=str, title=str, output=str, legend_position=str,
               conditions=str, show=begin.utils.tobool)
def main(*salmon_paths: "Path to salmon files",
         title: "The title of the graph"=
                "PCA of sample abundance over principal axis",
         output: "Path to the output file"="PCA.png",
         legend_position: "Position of the legend"="upper center",
         conditions: "Comma separated list of conditions per sample"=None,
         show: "Should the graph be displayed instead of saved ?"=False):
    """Plot a PCA on selected samples"""
    # Load datasets
    data = readSalmon(*salmon_paths)

    # Perform pca
    logging.debug("Computing components")
    nbc = len(salmon_paths)
    skpca = skd.PCA(n_components=nbc)

    # Prepare plots
    logging.debug("Fitting results")
    sktransform = skpca.fit_transform(data.T)
    results = pd.DataFrame(
        sktransform,
        columns=["PC_%i" % i for i in range(1, nbc+1, 1)],
        index=salmon_paths
    )

    logging.debug("Building graph")
    sns.set(style="darkgrid")
    try:
        results["Conditions"] = conditions.split(",")
        g = sns.FacetGrid(results, hue="Conditions")
    except AttributeError:
        g = sns.FacetGrid(results)

    g = g.map(plt.scatter, "PC_1", "PC_2")
    plt.title(title)
    try:
        frame = plt.legend(loc=legend_position, frameon=True).get_frame()
        frame.set_facecolor("white")
    except AttributeError:
        pass

    logging.debug("Plotting to %s" % ("standard output" if show else output))
    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
