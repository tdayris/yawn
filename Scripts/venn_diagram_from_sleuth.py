#!/usr/bin/python3
# conding: utf-8

import begin
import logging
import matplotlib.pyplot as plt
import pandas as pd
import traceback

from matplotlib_venn import venn2


@begin.convert(path=str, q_threashold=float, fc_threashold=float)
def read_curated(path, stat_threashold, fc_threashold):
    """
    Read a curated file:
    """
    logging.debug("Reading: %s" % path)
    data = pd.read_csv(path, sep="\t", header=0, index_col=0)
    logging.debug(data.head())
    fc_interest = [
        i for i in data.columns.tolist() if i.startswith("FC(Test_")
    ][0]
    data = data[["Q-Value", fc_interest]]
    data = data[data["Q-Value"] < stat_threashold]
    data = data[data[fc_interest] >= fc_threashold]
    logging.debug(data.head())
    logging.debug("Data shape is: %s" % str(data.shape))
    return set(data.index.tolist())


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(curated_sleuth_1=str, curated_sleuth_2=str, output=str,
               q_threashold=float, fc_threashold=float,
               show=begin.utils.tobool)
def main(curated_sleuth_1: "Path to a curated Sleuth file",
         set1_name: "Name of the set1",
         curated_sleuth_2: "Path to a second curated Sleuth file",
         set2_name: "Name of the second set",
         output: "Output prefix" = "Venn",
         q_threashold: "Threashold above which qval is significant" = 0.05,
         fc_threashold: "Threashold above which FC is significant" = 1,
         show: "Show graph instead of saving it to a file" = False):
    """
    This script builds a venn diagramm among two curated Sleuth files.

Example of curated Sleuth file:

GeneIdentifier	TranscriptIdentifier	Q-Value	stat_change	FC(Test_D1/Ref_D0)
RNA5S1	ENSG00000199352.1	1.6520570e-06	-2.987555871	-101.72037083
SPP1	ENSG00000118785.13	9.9955378e-05	1.9717884427	7.79610821734
EGR1	ENSG00000120738.7	9.9955325-05	2.3872604942	15.8468125203
CCR2	ENSG00000121807.5	9.9955325e-05	-2.683142349	-8.1815541162
    """
    sl1 = read_curated(curated_sleuth_1, q_threashold, fc_threashold)
    sl2 = read_curated(curated_sleuth_2, q_threashold, fc_threashold)

    venn_dict = {
        set1_name: sl1 - sl2,
        "%s & %s" % (set1_name, set2_name): sl1 & sl2,
        set2_name: sl2 - sl1
    }

    if len(venn_dict[set1_name]) == 0:
        set1_label = ""

    g = venn2(
        (len(sl1 - sl2), len(sl1 & sl2), len(sl2 - sl1)),
        set_labels=(
            set1_name if len(venn_dict[set1_name]) > 0 else "",
            set2_name if len(venn_dict[set2_name]) > 0 else ""
        )
    )

    g.get_patch_by_id("01").set_color("blue")
    g.get_patch_by_id("11").set_color("#DCDCDC")
    g.get_patch_by_id("10").set_color("red")

    if show:
        plt.show()
    else:
        plt.savefig("%s.png" % output, bbox_inches='tight')
        with open("%s.txt" % output, 'w') as outfile:
            outfile.write("%s: %s\n%s: %s\n%s: %s" % (
                set1_name,
                ", ".join(list(sl1 - sl2)),
                "%s & %s" % (set1_name, set2_name),
                ", ".join(list(sl1 & sl2)),
                set2_name,
                ", ".join(list(sl2 - sl1))
            ))
