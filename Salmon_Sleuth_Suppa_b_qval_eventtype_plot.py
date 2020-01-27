#!/usr/bin/python3.7
# conding: utf-8


import begin
import itertools
import logging
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels
import traceback

from collections import Counter
from matplotlib.colors import LogNorm
from matplotlib_venn import venn2
from pathlib import Path


@begin.convert(n1=str, n2=str, overwrite=begin.utils.tobool,
               show=begin.utils.tobool, stat_threshold=float,
               beta_threshold=float)
def combined_graphs(n1, n2, merged, overwrite=False, show=False,
                    beta_threshold=0.1, stat_threshold=0.01):
    pass


@begin.convert(n1=str, n2=str,
               overwrite=begin.utils.tobool,
               show=begin.utils.tobool,
               beta_threshold=float,
               stat_threshold=float)
def compare_frames(n1, n2, merged1, merged2, show=False, overwrite=False,
                   beta_threshold=0.1, stat_threshold=0.01):
    sdict = merged2.set_index("target_id").Significance.to_dict()
    merged1[f"{n2}_Status"] = [
        sdict[i] if i in sdict.keys() else "Below beta threshold"
        for i in merged1.target_id.tolist()
    ]
    logging.debug(Counter(merged1[f"{n2}_Status"]))

    outplot = Path(f"{n1}_vs_{n2}_Volcano.png")
    if (not outplot.exists()) or overwrite:
        sns.relplot(
            x="b",
            y="-Log10(Q-Value)",
            col="Event_Type",
            hue=f"{n2}_Status",
            hue_order=["Below significance threshold",
                       "Below beta threshold",
                       "Differentially Expressed"],
            col_order=['RI', 'A5', 'SE', 'A3', 'AL', 'AF', 'MX'],
            edgecolor="",
            palette="Blues",
            col_wrap=3,
            data=merged1
        )
        sns.despine(left=True)
        if show:
            plt.tight_layout()
            plt.show()
        else:
            plt.savefig(outplot, bbox_inches='tight')
        plt.clf()

    outplot = Path(f"{n1}_vs_{n2}_Volcano_significance_axis.png")
    if (not outplot.exists()) or overwrite:
        g = sns.lmplot(
            x="b",
            y="-Log10(Q-Value)",
            hue=f"{n2}_Status",
            hue_order=["Below significance threshold",
                       "Below beta threshold",
                       "Differentially Expressed"],
            palette="Blues",
            data=merged1,
            fit_reg=False,
            legend=False
        )
        g.map(
            plt.axhline,
            y=-np.log10(stat_threshold),
            ls=":",
            c=".5",
            markersize=0.1,
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
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = {
            l: h for l, h in zip(labels, handles)
            if not isinstance(h, matplotlib.lines.Line2D)
        }
        plt.legend(
            by_label.values(),
            by_label.keys(),
            bbox_to_anchor=(1.05, 1),
            loc=2, borderaxespad=0.
        )
        if show:
            plt.tight_layout()
            plt.show()
        else:
            plt.savefig(outplot, bbox_inches='tight')
        plt.clf()

    outplot = Path(f"{n1}_vs_{n2}_Venn.png")
    if (not outplot.exists()) or overwrite:
        de_targets1 = merged1["Significance"] == 'Differentially Expressed'
        de_targets1 = set(merged1[de_targets1].target_id)
        de_targets2 = merged2["Significance"] == 'Differentially Expressed'
        de_targets2 = set(merged2[de_targets2].target_id)
        venn2([de_targets1, de_targets2], [n1, n2])
        if show:
            plt.tight_layout()
            plt.show()
        else:
            plt.savefig(outplot, bbox_inches='tight')
        plt.cla()


@begin.convert(path=str, beta_threshold=float, stat_threshold=float)
def read_sleuth(path, beta_threshold=0.1, stat_threshold=0.01):
    columns_to_consider = ["qval", "b"]
    data = (pd.read_csv(path, sep="\t", index_col=False, header=0)
              .sort_values(by=["qval"], ascending=False)
              .dropna(how="all", subset=columns_to_consider))
    data["-Log10(Q-Value)"] = [
        0 if d == 0 else -np.log10(d) for d in data.qval
    ]
    data["Significance"] = list(classify(
        data.b, data.qval, beta_threshold, stat_threshold
    ))
    logging.debug(data.head())
    return data


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
@begin.convert(ioe_per_transcript=str, sleuth_paths=str,
               beta_threshold=float, stat_threshold=float,
               overwrite=begin.utils.tobool, show=begin.utils.tobool)
def main(ioe_per_transcript: "Path to a event per transcript files",
         *sleuth_paths: "Path to multiple Sleuth results",
         beta_threshold: "FC threshold" = 0.1,
         stat_threshold: "Significance threshold" = 0.01,
         overwrite: "Should outfile be overwritten or not?" = False,
         show: "Should graphs be plotted instead of saved" = False):
    """
    This script plots several graphs from Sleuth and Suppa outputs
    """
    sleuth_paths = [Path(i) for i in sleuth_paths]
    print(sleuth_paths)
    ioe_per_transcript = Path(ioe_per_transcript)

    ioe = pd.read_csv(
        str(ioe_per_transcript),
        sep="\t",
        index_col=False,
        header=0
    )
    logging.debug(ioe.head())

    sleuths = {
        i.parent: read_sleuth(str(i), beta_threshold, stat_threshold)
        for i in sleuth_paths
    }

    sns.set(style="darkgrid")
    for name, data in sleuths.items():
        outpath = Path(f"{name}_table.tsv")
        if (not outpath.exists()) or overwrite:
            merged = pd.merge(
                data,
                ioe,
                left_on="target_id",
                right_on="Transcript_ID",
                how="left"
            )
            logging.debug(merged.head())

            merged["GeneIdentifier"] = merged.ext_gene
            merged["stat_change"] = merged.b
            merged.set_index("target_id", inplace=True)
            merged = merged[
                merged["Significance"] == "Differentially Expressed"
            ]
            merged.to_csv(f"{name}_table.tsv", sep="\t")

    for n1, n2 in itertools.combinations(sleuths.keys(), 2):
        logging.debug(f"Working on: {n1} and {n2}")

        merged1 = pd.merge(
            sleuths[n1],
            ioe,
            left_on="target_id",
            right_on="Transcript_ID",
            how="left"
        )

        merged2 = pd.merge(
            sleuths[n2],
            ioe,
            left_on="target_id",
            right_on="Transcript_ID",
            how="left"
        )

        for name, merged in zip([n1, n2], [merged1, merged2]):
            outplot = Path(
                f"{name}_relplot_volcano_type_significance.png"
            )
            if (not outplot.exists()) or overwrite:
                sns.relplot(
                    x="b",
                    y="-Log10(Q-Value)",
                    col="Event_Type",
                    hue="Significance",
                    hue_order=["Below significance threshold",
                               "Below beta threshold",
                               "Differentially Expressed"],
                    col_order=['RI', 'A5', 'SE', 'A3', 'AL', 'AF', 'MX'],
                    edgecolor="",
                    palette="Blues",
                    # palette={
                    #     "Below significance threshold": "black",
                    #     "Differentially Expressed": "red",
                    #     "Below beta threshold": "grey"
                    # },
                    col_wrap=3,
                    data=merged
                )
                sns.despine(left=True)
                if show:
                    plt.tight_layout()
                    plt.show()
                else:
                    plt.savefig(outplot, bbox_inches='tight')
                plt.clf()

            outplot = Path(f"{name}_relplot_volcano_type_only.png")
            if (not outplot.exists()) or overwrite:
                sns.relplot(
                    x="b",
                    y="-Log10(Q-Value)",
                    hue="Event_Type",
                    hue_order=['RI', 'A5', 'SE', 'A3', 'AL', 'AF', 'MX'],
                    edgecolor="",
                    data=merged
                )
                sns.despine(left=True)
                if show:
                    plt.tight_layout()
                    plt.show()
                else:
                    plt.savefig(outplot, bbox_inches='tight')
                plt.clf()

        merged = pd.merge(
            sleuths[n1],
            sleuths[n2],
            on=["target_id", "ext_gene"],
            how="outer"
        )
        logging.debug(merged.head())

        combined_graphs(n1, n2, merged, overwrite, show,
                        beta_threshold, stat_threshold)

        compare_frames(n1, n2, merged1.copy(), merged2.copy(), show, overwrite,
                       beta_threshold, stat_threshold)
        compare_frames(n2, n1, merged2.copy(), merged1.copy(), show, overwrite,
                       beta_threshold, stat_threshold)
