#!/usr/bin/python3
# conding: utf-8

"""
Plot a pairwise scatterplot of multiple samples

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
import pandas as pd
import seaborn as sns
import traceback


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(salmon_paths=str, output=str, show=begin.utils.tobool)
def main(*salmon_paths: "Path to salmon files",
         output: "Path to the output file"="Pairwise_scatterplot.png",
         show: "Should the graph be displayed instead of saved ?"=False):
    """Plot a pairwise scatterplot of multiple samples"""
    # Load datasets
    frames = None
    for path in salmon_paths:
        logging.debug("Reading Salmon file: %s" % path)
        data = pd.read_csv(path, sep="\t", header=0, index_col=0)
        data = data[["TPM"]]
        data.columns = [path]
        try:
            frames = pd.merge(frames, data, left_index=True, right_index=True)
        except ValueError:
            frames = data

    # Build graph
    logging.debug("Building graph")
    sns.set(style="ticks", color_codes=True)
    g = sns.pairplot(frames,
                     diag_kind="kde",
                     diag_kws=dict(shade=True))

    logging.debug("Plotting to %s" % ("standard output" if show else output))
    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
