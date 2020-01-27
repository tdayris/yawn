#!/usr/bin/python3
# conding: utf-8

"""This file compares IOE file with a Sleuth file"""

import argparse
import begin
import logging as lg
import os.path as op
import pandas as pd
import pygal as pg
import traceback


@begin.convert(path=str)
def readSleuth(path: "(str)", qval: "(float)"=0.001, diff_exp: "(bool)"=False):
    """
    Return a pandas.DataFrame from a sleuth file.
    Remove equelly expressed transcripts if requested.

    Parameters:
        path     (str)  : Path to a Sleuth tsv file
        qval     (float): p-valu/q-value used to select diff. exp. transcripts
        diff_exp (bool) : weather to trucate or not the DataFrame to
                          differentially expressed transcripts only,
                          using the given q-value

    Return:
        pandas.DataFrame
    """
    lg.debug("Reading Sleuth-TSV file: %s" % path)
    data = pd.read_csv(path, sep="\t", header=0, index_col=0)
    # Sleuth can have string values among its result table
    # Thiat's why errors are ignored in the next command.
    # String values won't be changed, while nummeric
    # will be correctly casted.
    data = data.apply(pd.to_numeric, errors="ignore")

    if diff_exp:
        lg.debug("Data are reduced to diff. exp. only: qval <= %f" % qval)
        data = data[data["qval"] <= qval]

    lg.debug(data.head())
    return data


@begin.convert(path=str, gene=str)
def readIOE(path: "(str)", gene: "(str)"=""):
    """
    Return a pandas.DataFrame from an ioe file.
    Keeps only the requested gene if requested.

    Parameters:
        path (str): Path to an IOE tsv file
        gene (str): The gene ID

    Return:
        pandas.DataFrame
        """
    lg.debug("Reading IOE-TSV file: %s" % path)
    # IOE's 3rd column is the identifier of the junction. It is unique
    # and contains many information repreated along the rest of the columns.
    data = pd.read_csv(path, sep="\t", header=0, index_col=2)
    lg.debug(data.head())
    return data


@begin.convert(path=str)
def readDPSI(path: "(str)"):
    """
    Return a pandas.DataFrame from a DPSI file.

    Parameters:
        path (str): Path to a DPSI tsv file

    Return:
        pandas.DataFrame
    """
    lg.debug("Reading DPSI-TSV file: %s" % path)
    data = pd.read_csv(path, sep="\t", header=0, index_col=0)
    # Every value in this table should be int or float (weather Nan or not)
    # If any value does not follow this sheme, then TypeError will be raised.
    data = data.apply(pd.to_numeric, errors="raise")
    lg.debug(data.head())
    return data


def getTrFromIOE(ioe: "(pandas.DataFrame)",
                 column: "(str)"="alternative_transcripts",
                 drop_duplicate: "(bool)"=True):
    """
    Return a pandas.Series of transcripts identifiers from an IOE file

    """
    pass


@begin.convert(drop_duplicate=begin.utils.tobool)
def getGeneFromJunctions(junctions: "(pandas.Series)",
                         drop_duplicate: "(bool)"=True):
    """
    Return a pandas.Serie of gene identifiers related to the given junctions

    Parameters:
        junctions      (pandas.Series): A Serie of junctions' names
        drop_duplicate (bool)         : Weather or not to frop duplicates

    Return:
        pandas.Series"""
    # Split the column according to commas, format-it into a 1d Serie
    genes = junctions.str.split(";").apply(pd.Series, 1)

    # Genes names are in column number 0. The rest of the dataframe is
    # composed of irrelevant junction information
    genes = genes[genes.columns.tolist()[0]]

    if drop_duplicate:
        # Inplace dropping to inscrease speed
        genes.drop_duplicates(keep="first", inplace=True)

    lg.debug(genes.head())
    return genes


@begin.start(formatter_class=argparse.RawTextHelpFormatter, cmd_delim='--')
@begin.logging
@begin.tracebacks
@begin.convert(sleuth_path=str, ioe_path=str, colname=str, analyses=str,
               output=str, on_genes=begin.utils.tobool,
               on_transcripts=begin.utils.tobool, qval=float,
               drop_equally_expressed=begin.utils.tobool)
def compareSleuthIOE(sleuth_path: "(str)", ioe_path: "(str)",
                     output: "(str)", analyses: "(str)"="all",
                     colname: "(str)"="alternative_transcripts",
                     qval: "(float)"=0.001, on_genes: "(bool)"=False,
                     on_transcripts: "(bool)"=True,
                     drop_equally_expressed: "(bool)"=True):
    """Compare both Sleuth and IOE tsv-file.

    You may analyse both features : Genes and/or transcripts.
    The following alanyses are possible:
        Conjunction            AND gate    conjunction
        Disjunction            OR gate     disjunction
        Exclusive Disjunction  XOR gate    exclusive_disjunction
        Negation               NAND gate   negation
        All of the above                   all

    Parameters:
        sleuth_path            (str):  A sleuth tsv-file path
        ioe_path               (str):  An IOE tsv-file path
        output                 (str):  The stem name of the output file
        colname                (str):  The name of the column containing
                                       genes id in <ioe>
        analyses               (str):  The name of the analyse that is to be
                                       performed qval (float): The q-value
                                       threashold to identify diff. exp.
                                       features
        on_genes               (bool): Weather or not to perform analyses on
                                       genes, can be used with <on_transcripts>
        on_transcripts         (bool): Weather ot not to perform analyses
                                       on transcripts (can be used with
                                       <on_genes>)
        drop_equally_expressed (bool): Weather to drop or not equally expressed
                                       features

    Return:
        None.
        Writes TSV files on the disk.
    """
    available_analyses = {
        # By using this dict, strings predicates are replaced
        # by bitwise int comparison. It save one if/else block,
        # and is quicker to answer.
        "conjunction": 1,
        "disjunction": 2,
        "exclusive_disjunction": 4,
        "negation": 8,
        "all": 15,
    }

    try:
        analyses = available_analyses[analyses]
        lg.debug("Analyse recognized.")
    except KeyError:
        msg = ("Unavailable analyse: %s\n" % analyses,
               "Expected any of: %s" % list(available_analyses.keys()))
        lg.error(msg)
        raise KeyError(msg)

    # Loading data tables into memory
    sleuth = readSleuth(sleuth_path, qval, drop_equally_expressed)
    ioe = readIOE(ioe_path)

    # Getting the transcripts column, splitting according to commas,
    # then reshaping it into a pandas series, and stacking the whole
    # thing into a 1d Series and keeping the original event name as
    # a meta-data.
    tr_id = (set(ioe[colname].str.split(",").apply(pd.Series, 1).stack())
             if on_transcripts else None)

    gene_id = (set(getGeneFromJunctions(pd.Series(ioe.index)))
               if on_genes else None)

    lg.debug("Subsets extracted")

    name, ext = op.splitext(output)
    if analyses & 1:  # Conjunction or all
        # tr_id and gene_id are either sets or NoneType. The following
        # will rase an error if and only if tr_id and/or gene_id are NoneType.
        # It's better to ask forgiveness than permission, in python !
        try:
            # Retrieve information
            data = sleuth[
                sleuth.index.isin(set(sleuth.index.tolist()) & tr_id)
            ]

            # Save it
            data.to_csv(
                "%s_tr_conjunction%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Transcript Conjunction: %s" % str(data.shape))
        except TypeError:
            lg.debug("No Transcript Conjunction requested.")

        try:
            # Retrieve information
            data = sleuth[
                sleuth["ens_gene"].isin(
                    set(sleuth["ens_gene"].tolist()) & gene_id
                )
            ]

            # Save it
            data.to_csv(
                "%s_gene_conjunction%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Gene Conjunction: %s" % str(data.shape))
        except TypeError:
            lg.debug("No Gene Conjunction requested.")

    if analyses & 2:  # Disjunction or all
        try:
            # Retrieve information
            data = sleuth[
                sleuth.index.isin(set(sleuth.index.tolist()) | tr_id)
            ]

            # Save it
            data.to_csv(
                "%s_tr_disjunction%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Transcript Disjunction: %s" % str(data.shape))
        except TypeError:
            lg.debug("No Transcript Disjunction requested.")

        try:
            # Retrieve information
            data = sleuth[
                sleuth["ens_gene"].isin(
                    set(sleuth["ens_gene"].tolist()) | gene_id
                )
            ]

            # Save it
            data.to_csv(
                "%s_gene_disjunction%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Gene Disjunction: %s" % str(data.shape))
        except TypeError:
            lg.debug("No Gene Disjunction requested.")

    if analyses & 4:  # Exclusive Disjunction or all
        try:
            # Retrieve information
            data = sleuth[
                sleuth.index.isin(set(sleuth.index.tolist()) ^ tr_id)
            ]

            # Save it
            data.to_csv(
                "%s_tr_exlusive_disjunction%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Transcript Exclusive Disjunction"
                     ": %s" % str(data.shape))
        except TypeError:
            lg.debug("No Transcript Exclusive Disjunction requested.")

        try:
            # Retrieve information
            data = sleuth[
                sleuth["ens_gene"].isin(
                    set(sleuth["ens_gene"].tolist()) ^ gene_id
                )
            ]

            # Save it
            data.to_csv(
                "%s_gene_exlusive_disjunction%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Gene Exclusive Disjunction Saved %s" % str(data.shape))
        except TypeError:
            lg.debug("No Gene Exclusive Disjunction requested.")

    if analyses & 8:  # Negation or all
        try:
            # Retrieve information
            data = sleuth[
                sleuth.index.isin(set(sleuth.index.tolist()) - tr_id)
            ]

            # Save it
            data.to_csv(
                "%s_tr_negation%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Transcript Negation Saved %s" % str(data.shape))
        except TypeError:
            lg.debug("No Transcript Negation requested.")

        try:
            # Retrieve information
            data = sleuth[
                sleuth["ens_gene"].isin(
                    set(sleuth["ens_gene"].tolist()) - gene_id
                )
            ]

            # Save it
            data.to_csv(
                "%s_gene_negation%s" % (name, ext),
                sep="\t"
            )

            # Return information if verbose
            lg.debug("Gene Negation Saved %s" % str(data.shape))
        except TypeError:
            lg.debug("No Gene Negation requested.")
