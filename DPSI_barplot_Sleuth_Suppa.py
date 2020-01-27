#!/usr/bin/python3
# conding: utf-8

"""
Plot stacked bar of DPSI ASE usage across multiple sample

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
import pandas as pd
import seaborn as sns
import traceback


@begin.convert(paths=str)
def readDPSI(*paths: "Path to dpsi files"):
    """Return a DataFrame from a dpsi file"""
    frames = []
    for path in paths:
        logging.debug("Reading DPSI file: %s" % path)

        # Load dataframe
        data = pd.read_csv(path, sep="\t", header=0, index_col=None)

        # Split columns
        try:
            cols = (pd.DataFrame(data["Position"].str
                                                 .split(";")
                                                 .apply(pd.Series)))
            del data["Position"]
        except KeyError:
            cols = (pd.DataFrame(data["Positions"].str
                                                  .split(";")
                                                  .apply(pd.Series)))
            del data["Positions"]

        cols.columns = ["Gene", "Position"]
        data.columns = ["DPSI", "pval"]
        data = data.join(cols)

        cols = pd.DataFrame(data["Position"].str.split(":").apply(pd.Series))
        try:
            cols.columns = ["ASE Type", "Sequence",
                            "Event1", "Event2", "Strand"]
        except ValueError:
            try:
                cols.columns = ["ASE Type", "Sequence", "Event1_Start",
                                "Event1_Stop", "Event2_Start", "Event2_Stop",
                                "Strand"]
            except ValueError:
                cols.columns = ["ASE Type", "Sequence", "Start",
                                "Intron", "Stop", "Strand"]
        data = data.join(cols)

        # Define index col
        data.set_index("Position", inplace=True)

        frames.append(data)
    return pd.concat(frames)


@begin.convert(path=str, pvalue=float, qvalue=float,
               is_gene=begin.utils.tobool)
def readSleuth(path: "Path to a Sleuth file",
               pvalue: "P-value cutoff"=0.01,
               qvalue: "Q-Value cutoff"=0.05,
               is_gene: "Is sleuth a DGE exp ?"=False):
    """Return a DataFrame from a Sleuth file"""
    logging.debug("Reading Sleuth file: %s" % path)

    # Load DataFrame
    data = pd.read_csv(path, sep="\t", index_col=0, header=0)

    # Reduce dataframe
    cols = ("target_id" if is_gene else "ens_gene")
    data = data[data["pval"] <= pvalue]
    data = data[data["qval"] <= qvalue]
    data = data[cols].drop_duplicates(inplace=True)

    return data


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(dpsi_path=str, sleuth_path=str, title=str, output=str,
               pvalue=float, qvalue=float, is_gene=begin.utils.tobool,
               show=begin.utils.tobool)
def main(sleuth_path: "Path to a Sleuth file",
         *dpsi_path: "Path to multiple dpsi files",
         title: "Title of the graph"="Variation of the inclusion of "
                                     "junctions across samples",
         output: "Path to the output file"="DPSI_ASE.png",
         pvalue: "P-value cutoff"=0.01,
         qvalue: "Q-Value cutoff"=0.05,
         is_gene: "Is sleuth a DGE exp ?"=False,
         show: "Should the graph be displayed instead of saved ?"=False):
    """Create a Barplot from multiple dpsi files, considering Sleuth DTE"""
    genes_index = readSleuth(sleuth_path, pvalue, qvalue, is_gene)
    data = readDPSI(*dpsi_path)
    logging.debug(data.head())

    sns.set(style="whitegrid")
    ax = sns.violinplot(x="ASE Type",
                        y="DPSI",
                        data=data,
                        palette="Set3",
                        linewidth=1,
                        bw=.2)
    ax.set_title(title)

    sns.despine(left=True, bottom=True)
    if show:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')
