#!/usr/bin/python3


import begin
import collections
import logging
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import pandas as pd
import re
import seaborn as sns
import sys
import traceback

sys.path.append("/home/tdayris/Documents/Developpement/STRonGR/scripts/")

import target_clustered_heatmap as tch
import Clustermap as clmap


@begin.convert(path=str, gene=begin.utils.tobool)
def quant_reader(path: "Path to a tsv file",
                 gene: "Is it a gene quantification?" = False) \
        -> pd.DataFrame:
    """
    The quantification specific tsv file reader
    """
    id_name = "Gene_ID" if gene is True else "Transcript_ID"
    data = pd.read_csv(path, sep="\t", index_col=0, header=0, na_values="")
    # logging.debug(data.head())
    return data


@begin.convert(path=str)
def sleuth_reader(path: "Path to a tsv file",
                  condition: "The condition name" = "",
                  gene: "Is it a DGE?" = False) -> pd.DataFrame:
    """
    The sleuth specific tsv file reader
    """
    id_name = "Gene_ID" if gene is True else "Transcript_ID"
    data = pd.read_csv(
        path,
        sep="\t",
        index_col=([11, 0] if gene is True else [12, 11, 0]),
        header=0,
        dtype={
            0: str,
            1: pd.np.float,
            2: pd.np.float,
            3: pd.np.float,
            4: pd.np.float,
            5: pd.np.float,
            6: pd.np.float,
            7: pd.np.float,
            8: pd.np.float,
            9: pd.np.float,
            10: pd.np.float,
            11: str,
            12: str
        },
        na_values=""
    )
    if gene is True:
        data = data.reset_index()[["ext_gene", "target_id", "qval", "b"]]
        data.columns = ["Gene_Name", "Gene_ID",
                        "QValue_" + condition, "BetaValue_" + condition]
        data.set_index(["Gene_Name", "Gene_ID"], inplace=True)
    else:
        data = data.reset_index()[[
            "ext_gene", "ens_gene", "target_id", "qval", "b"
        ]]
        data.columns = [
            "Gene_Name", "Gene_ID", "Transcript_ID",
            "QValue_" + condition, "BetaValue_" + condition
        ]
        data.set_index(["Gene_Name", "Gene_ID", "Transcript_ID"], inplace=True)

    # logging.debug(data.head())
    return data


@begin.convert(path=str)
def read_t2g(path: "Path to a T2G tsv file",
             genes: "Return genes only") -> pd.DataFrame:
    """
    T2G are three columned tsv files with:

    1: ens_gene
    2: ens_tr
    3: ext_gene
    """
    t2g = pd.read_csv(path, sep="\t", header=None, index_col=1)
    t2g.columns = ["Ensembl_ID", "Hugo_ID"]
    if genes:
        t2g = (t2g.reset_index()[["Ensembl_ID", "Hugo_ID"]]
                  .drop_duplicates()
                  .set_index("Ensembl_ID"))
    # logging.debug(t2g.head())
    return t2g


def variation_level(diffexp: "Pandas DataFrame build within preprocessing",
                    condition: "Condition name",
                    threshold: "The stat value threshold"):
    """
    Generates the variation level column
    """
    zipped_cols = zip(
        diffexp["QValue_%s_DTE" % condition],
        diffexp["BetaValue_%s_DTE" % condition],
        diffexp["QValue_%s_DGE" % condition],
        diffexp["BetaValue_%s_DGE" % condition]
    )
    for qtr, btr, qge, bge in zipped_cols:
        value = []
        if float(qtr) <= float(threshold):
            if float(btr) < 0:
                value.append("Up")
            else:
                value.append("Down")
        else:
            value.append("EE")
        if float(qge) <= float(threshold):
            if float(bge) < 0:
                value.append("Up")
            else:
                value.append("Down")
        else:
            value.append("EE")
        yield ":".join(value)


def infer_splice(diffexp: "Pandas DataFrame build within preprocessing",
                 colname: "Column of interest"):
    """
    Generates the splicing inferance
    """
    for v in diffexp[colname]:
        v = v.split(":")
        if v[0] == v[1]:
            yield False
        else:
            yield True


@begin.convert(diff_file=str, quant_file=str, output_prefix=str,
               qval_threshold=float)
def clustermap(diff_file: "Path to the DE file",
               quant_file: "Path to the quant file",
               output_prefix: "Prefix of the output file",
               condition_array: "The condition array",
               qval_threshold: "The qvalue threshold"):
    """
    Saves clstered heatmaps
    """
    tch.main(
        diff_file,
        quant_file,
        *condition_array,
        qval=qval_threshold,
        max_target=100,
        output=output_prefix + "_TargetClusteredHeatmap.png",
        blue=False,
        show=False
    )

    clmap.main(
        quant_file,
        *condition_array,
        title="Clustered Heatmap",
        output=output_prefix + "_SampleClusteredHeatmap.png",
        xlabel_rotation=0,
        ylabel_rotation=90,
        row_cluster=True,
        col_cluster=True,
        row_condition_color=True,
        col_condition_color=True,
        show=False,
        square=False,
        robust=False
    )


@begin.convert(path=str)
def read_psi(path: "Path to PSI file",
             condition_name: "Name of the current condition"):
    """Read a PSI file and return a DataFrame"""

    # Load DataFrame
    data = pd.read_csv(path, sep="\t", header=0, index_col=0)
    data.reset_index(inplace=True)
    data.columns = ["Index"] + \
                   ["PSI_%s_%s" % (condition_name, i)
                    for i in data.columns.tolist()[1:]]

    # Split first column
    cols = pd.DataFrame(data["Index"].str.split(";", 1).tolist())
    del data["Index"]

    cols.columns = ["Gene_ID", "Transcript_ID"]
    data = data.join(cols)
    # logging.debug(data.head())

    return data.set_index(["Gene_ID", "Transcript_ID"])


@begin.convert(transcript_quant=str, gene_quant=str, transcripts_diff=str,
               genes_diff=str, qval_threshold=str)
def run(transcript_quant: "Path to a normalized quantification table",
        gene_quant: "Path to a normalized gene quantification table",
        transcripts_diff: "Path to transcript level differential analysis",
        genes_diff: "Path to gene level differential analysis",
        transcript_psi: "Path to PSI file from Suppa",
        sample_condition: "Mapping Sample/Condition",
        qval_threshold: "The Q-Value threshold for significance" = 0.05,
        condition_name: "The condition name" = "",
        cluster_splice=True,
        plot_isoform_psi=True) \
            -> pd.DataFrame:
    """
    Create variation level column
    """
    # Pre-processing
    logging.debug("Loading Sleuth results")
    dge = sleuth_reader(genes_diff, condition_name, True)
    dte = sleuth_reader(transcripts_diff, condition_name, False)

    diffexp = pd.merge(
        dte.reset_index(),
        dge.reset_index(),
        on=["Gene_Name", "Gene_ID"],
        how="left",
        suffixes=["_DTE", "_DGE"]
    ).set_index(["Gene_Name", "Gene_ID", "Transcript_ID"])
    logging.debug("Original shape: %s, %s" % (str(dte.shape), str(dge.shape)))
    logging.debug("Merged shape: %s" % str(diffexp.shape))

    vlevel = condition_name + "_Variation_Level"
    splice = condition_name + "_Splicing"

    diffexp[vlevel] = list(
        variation_level(diffexp, condition_name, qval_threshold)
    )
    logging.debug("Regulation level overview: %s" % str(
        collections.Counter(diffexp[vlevel]))
    )

    diffexp[splice] = list(
        infer_splice(diffexp, vlevel)
    )
    logging.debug("Alternative Splicing: %s" % str(
        collections.Counter(diffexp[splice]))
    )

    splice_list = diffexp[diffexp[splice]].reset_index()

    logging.debug("Loading Quant Tr results")
    quant_tr = quant_reader(transcript_quant)
    cond_sample = {
        k: v
        for k, v in zip(set(quant_tr.columns), sample_condition.split(","))
    }
    quant_tr.columns = [
        "%s_Quant_Tr_%s" % (condition_name, i)
        for i in quant_tr.columns.tolist()
    ]
    # for cond in set(sample_condition.split(",")):
        # cond_cols = [i for i, j in cond_sample.items() if j == cond]
        # quant_tr["Sum_Tr_" + cond] = quant_tr.loc[:, quant_tr.columns.isin(cond_cols)].sum(1)

    logging.debug("Transcript quant shape: %s" % str(quant_tr.shape))
    splice_tr = quant_tr[quant_tr.index.isin(splice_list["Transcript_ID"])]
    logging.debug("Transcript splice shape: %s" % str(splice_tr.shape))
    quant_tr.to_csv(
        condition_name + "_Quantification_Spliced_Tr.tsv",
        sep="\t"
    )
    condition_array = [condition_name]*len(quant_tr.columns)
    if cluster_splice is True:
        logging.debug("Making clustermap out of quant_tr")
        clustermap(
            transcripts_diff,
            condition_name + "_Quantification_Spliced_Tr.tsv",
            condition_name + "_Tr",
            [condition_name]*len(quant_tr.columns),
            qval_threshold
        )

    logging.debug("Loading Quant Ge results")
    quant_ge = quant_reader(gene_quant)
    quant_ge.columns = [
        "%s_Quant_Ge_%s" % (condition_name, i)
        for i in quant_ge.columns.tolist()
    ]
    cond_sample = {
        k: v
        for k, v in zip(set(quant_ge.columns), sample_condition.split(","))
    }
    # for cond in set(sample_condition.split(",")):
        # quant_ge["Sum_Ge_" + cond] = quant_ge[
        #     [i for i, j in cond_sample.items() if j == cond]
        # ].sum(1)
    logging.debug("Gene quant shape: %s" % str(quant_ge.shape))
    splice_ge = quant_ge[quant_ge.index.isin(splice_list["Gene_ID"])]
    logging.debug("Gene splice shape: %s" % str(splice_ge.shape))
    quant_ge.to_csv(
        condition_name + "_Quantification_Spliced_Ge.tsv",
        sep="\t"
    )
    condition_array = [condition_name]*len(quant_ge.columns)
    if cluster_splice is True:
        logging.debug("Making clustermap out of quant_tr")
        clustermap(
            genes_diff,
            condition_name + "_Quantification_Spliced_Ge.tsv",
            condition_name + "_GE",
            [condition_name]*len(quant_ge.columns),
            qval_threshold
        )
        logging.debug("Done")

    diffexp = pd.merge(
        diffexp.reset_index(),
        quant_tr.reset_index(),
        left_on=["Transcript_ID"],
        right_on=["index"],
        how='left'
    ).set_index(["Gene_Name", "Gene_ID", "Transcript_ID"])
    del diffexp["index"]

    diffexp = pd.merge(
        diffexp.reset_index(),
        quant_ge.reset_index(),
        left_on=["Gene_ID"],
        right_on=["index"],
        how='left'
    ).set_index(["Gene_Name", "Gene_ID", "Transcript_ID"])
    del diffexp["index"]

    psi = read_psi(transcript_psi, condition_name)
    diffexp = pd.merge(
        diffexp.reset_index(),
        psi.reset_index(),
        on=["Gene_ID", "Transcript_ID"],
        how='left'
    ).set_index(["Gene_Name", "Gene_ID", "Transcript_ID"])

    psi = (psi.reset_index()
              .melt(id_vars=["Gene_ID", "Transcript_ID"],
                    var_name=["Sample_ID"],
                    value_name="PSI"))
    cond_sample = {
        k: v
        for k, v in zip(set(psi["Sample_ID"]), sample_condition.split(","))
    }

    psi["Condition"] = [cond_sample[k] for k in psi["Sample_ID"]]
    psi = pd.merge(
        psi,
        diffexp.reset_index()[["Transcript_ID", splice, vlevel]],
        on="Transcript_ID",
        how="left"
    )
    psi.set_index(["Gene_ID", "Transcript_ID"], inplace=True)
    printed_list = []
    if plot_isoform_psi is True:
        sns.set(style="darkgrid")
        for gene in set(psi[psi[splice]].index.get_level_values("Gene_ID")):
            subset = (psi.loc[psi.index
                                 .get_level_values("Gene_ID")
                                 .str
                                 .startswith(gene)]
                         .dropna(how='any'))
            if len(subset) / len(sample_condition.split(",")) <= 2:
                logging.debug("Skipping %s: too view transcripts" % gene)
                continue

            if all(all(j <= 0.3 for j in pd.Series(i[1].std()))
                   for i in subset.groupby("Transcript_ID")):
                logging.debug("Skipping %s: too tight variations" % gene)
                continue

            logging.debug("Working on %s" % gene)

            printed_list.append(gene)

            g = sns.PairGrid(
                subset.reset_index(),
                y_vars=["Transcript_ID"],
                x_vars=["PSI"],
                height=10,
                aspect=.25,
                hue="Condition"
            )
            g.map(
                sns.stripplot,
                size=10,
                orient="h",
                linewidth=1,
                jitter=False
            )
            g.add_legend()
            g.savefig(condition_name + "_" + gene + "_DotPlot_PSI_Isoforms.png",
                      bbox_inches='tight')
            logging.debug("Saved: %s" % gene + "_DotPlot_PSI_Isoforms.png")

    if not printed_list == []:
        print(" ".join(printed_list))
    psi.to_csv("PSI_melter_" + condition_name + ".tsv", sep="\t")

    # logging.debug(diffexp.head())
    return diffexp


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(plot_isoform_psi=begin.utils.tobool,
               cluster_splice=begin.utils.tobool)
def main(transcript_quant: "Path to normalized quantification "
                           "tables (comma separated)",
         gene_quant: "Path to a normalized gene quantification "
                     "tables (comma separated)",
         transcripts_diff: "Path to transcript level differential "
                           "analysis tables (comma separated)",
         genes_diff: "Path to gene level differential analysis "
                     "tables (comma separated)",
         condition_names: "Comma separated list of conditions",
         transcript_psi: "Path to transcript psi files",
         appris_data: "Path to appris annotation file",
         *sample_condition: "Multiple comma separated list"
                            " of condition per sample",
         plot_isoform_psi: "Perform a Paired Stripplot of Isoforms" = True,
         cluster_splice: "Perform a Clustermap of protentially "
                         "spliced isoforms" = True,
         # tr2gene: "Path to a transcript to gene tsv table",
         qval_threshold: "The Q-Value threshold for significance" = 0.05) \
            -> None:
    """
    This script performs classical table merging for DTE/DGE integration with:

    1- Salmon
    2- Sleuth
    3- Python
    """
    sns.set(style="darkgrid")

    condition_zip = zip(
        transcript_quant.split(","),
        gene_quant.split(","),
        transcripts_diff.split(","),
        genes_diff.split(","),
        condition_names.split(","),
        transcript_psi.split(","),
        sample_condition
    )
    complete_table = None
    for i, j, k, l, m, n, o in condition_zip:
        data = run(i, j, k, l, n, o, qval_threshold,
                   m, cluster_splice, plot_isoform_psi)
        try:
            tmp = pd.merge(
                complete_table,
                data,
                left_index=True,
                right_index=True,
                how="right"
            )
            logging.debug(
                "Original shape: %s, %s" % (
                    str(complete_table.shape),
                    str(data.shape))
                )
            logging.debug("Merged shape: %s" % str(tmp.shape))
        except ValueError:
            complete_table = data
        else:
            complete_table = tmp
            del tmp
        finally:
            del data

    complete_table.to_csv(
        "Complete_table_%s.tsv" % "_".join(condition_names.split(',')),
        sep="\t"
    )

    complete_table["Web_Summary"] = [
        "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s" % idx[1]
        for idx in complete_table.index
    ]

    complete_table.reset_index(inplace=True)
    complete_table["Transcript_ID"] = [
        re.sub('\..*$', '', i)
        for i in complete_table["Transcript_ID"]
    ]
    complete_table["Gene_ID"] = [
        re.sub('\..*$', '', i)
        for i in complete_table["Gene_ID"]
    ]
    complete_table.set_index(["Gene_Name", "Gene_ID", "Transcript_ID"])

    appris = pd.read_csv(
        appris_data,
        sep="\t",
        header=None,
        index_col=[0, 1]
    )
    appris.columns = ["Protein_ID", "APPRIS"]

    complete_table = pd.merge(
        complete_table.reset_index(),
        appris.reset_index(),
        left_on=["Gene_ID", "Transcript_ID"],
        right_on=[0, 1],
        how='left'
    )
    del complete_table[0]
    del complete_table[1]
    del complete_table["index"]
    complete_table.set_index(["Gene_Name", "Gene_ID", "Transcript_ID"])
    appris_values = ["PRINCIPAL:1", "PRINCIPAL:2", "PRINCIPAL:3",
                     "PRINCIPAL:4", "PRINCIPAL:5", "ALTERNATIVE:1",
                     "ALTERNATIVE:2", "MINOR"]
    complete_table["APPRIS"] = [
        value if value in appris_values else "-"
        for value in complete_table["APPRIS"]
    ]

    complete_table.to_csv(
        "Complete_table_%s.tsv" % "_".join(condition_names.split(',')),
        sep="\t"
    )

    tmp = complete_table = complete_table[
        [i for i in complete_table.columns
         if i.endswith("_Splicing") or
         (i.startswith("QValue_") and i.endswith("_DTE"))]
    ]

    for cond in condition_names.split(","):
        splice = cond + "_Splicing"
        subset = subset[subset[splice]]
        subset = subset[[i for i in subset.columns if "_%s_" % cond in i]]
        # logging.debug(subset.head())
        subset = subset.reset_index()[["Gene_ID", "QValue_%s_DTE" % cond]]
        subset.columns = ["ext_gene", "qval"]
        subset.to_csv("%s_upsetR.tsv")

    # Now, for Transcript then for Genes
    # Make hist of non EE target per condition
    # Make venn of non EE target per condition

    # Load Suppa data
    # Extract targets event type
    # Make Hist of event types count
