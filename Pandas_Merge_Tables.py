#!/usr/bin/python3
# conding: utf-8

"""
Join multiple Sleuth files, rename columns and stack them on demand

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
import pandas as pd
import traceback

from pathlib import Path


@begin.convert(salmon=begin.utils.tobool, sleuth=begin.utils.tobool)
def swith(salmon, sleuth) -> pd.DataFrame:
    """
    Calls the right function
    """
    if salmon:
        return salmon_reader
    elif sleuth:
        return sleuth_reader
    return generic_reader


@begin.convert(path=str)
def read_t2g(path: "Path to a T2G tsv file") -> pd.DataFrame:
    """
    T2G are three columned tsv files with:

    1: ens_gene
    2: ens_tr
    3: ext_gene
    """
    t2g = pd.read_csv(path, sep="\t", header=None, index_col=1)
    t2g.columns = ["Ensembl_ID", "Hugo_ID"]
    return t2g


@begin.convert(path=str)
def generic_reader(path: "Path to a tsv file") -> pd.DataFrame:
    """
    The most generic tsv file reader
    """
    return pd.read_csv(path, sep="\t", header=0, index_col=0)


@begin.convert(path=str)
def salmon_reader(path: "Path to a tsv file") -> pd.DataFrame:
    """
    The salmon specific tsv file reader
    """
    return pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        dtype={
            0: str,
            1: pd.np.float,
            2: pd.np.float,
            3: pd.np.float,
            4: pd.np.float
        },
        na_values=""
    )


@begin.convert(path=str)
def sleuth_reader(path: "Path to a tsv file") -> pd.DataFrame:
    """
    The sleuth specific tsv file reader
    """
    return pd.read_csv(
        path,
        sep="\t",
        index_col=[12, 11, 0],
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


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(path=str, col=str, output=str, prefix=str, suffix=str,
               index_label=begin.utils.tobool, stack=begin.utils.tobool,
               salmon=begin.utils.tobool, sleuth=begin.utils.tobool)
def main(*paths: "Multiple paths to TSV files",
         col: "Name of the column that it to be joined"="TPM",
         output: "Path to the output file"="Frames.tsv",
         prefix: "The prefix to be removed to get "
                 "sample name from its file name"="",
         suffix: "The suffix to be removed to get "
                 "sample name from its file name"="",
         index_label: "Should the index be labelled ?"=True,
         stack: "Should dataset be stacked ?"=False,
         salmon: "All files are Salmon files"=False,
         sleuth: "All files are Sleuth files"=False,
         drop_null: "Remove unused targets"=False,
         drop_na: "Remove full NaN lines"=False,
         tr2g: "Path to a Transcript to Gene table"=None) -> None:
    """
    This script reads several tsv files one after another. It relies over the
    follwing assumptions:

    1. The lines starting with '#' are comments.
    2. The first line without '#' are headers.
    3. The first column contains unique keys.

    The options:
    --salmon, --sleuth, do not take these previous assumptions into account.
    """
    merged_frame = None
    destin = Path(output)
    parser = swith(salmon, sleuth)
    if destin.exists():
        raise FileExistsError("Output file already exists: %s" % str(destin))

    for path in paths:
        logging.debug("Working on %s" % path)
        source = Path(path)
        if not source.exists():
            raise FileNotFoundError("Could not find: %s" % str(source))

        data = swith(salmon, sleuth)(str(source))
        sample_id = str(source)
        sample_id = sample_id if suffix == "" else sample_id[:-len(suffix)]
        sample_id = sample_id if prefix == "" else sample_id[len(prefix):]
        logging.debug("%s's shape: %s" % (sample_id, str(data.shape)))

        data = data[[col]]
        data.columns = [sample_id]

        try:
            if sleuth:
                merged_frame = pd.merge(
                    merged_frame.reset_index(),
                    data.reset_index(),
                    on=["ext_gene", "ens_gene", "target_id"]
                ).set_index(["ext_gene", "ens_gene", "target_id"])
            else:
                merged_frame = pd.merge(
                    merged_frame,
                    data,
                    left_index=True,
                    right_index=True
                )
        except ValueError:
            merged_frame = data
        except AttributeError:
            merged_frame = data

    logging.debug("Frames merged.")
    logging.debug(merged_frame.shape)
    logging.debug(merged_frame.head())

    if drop_null:
        logging.debug("Removing null lines")
        logging.debug("Shape: %s" % str(merged_frame.shape))
        merged_frame = merged_frame.loc[~(merged_frame == 0).all(axis=1)]
        logging.debug("Shape: %s" % str(merged_frame.shape))

    if drop_na:
        logging.debug("Removing NaN lines")
        logging.debug("Shape: %s" % str(merged_frame.shape))
        merged_frame.dropna(axis=0, how="all", inplace=True)
        logging.debug("Shape: %s" % str(merged_frame.shape))

    if stack:
        logging.debug("Stacking merged_frame")
        merged_frame = pd.DataFrame(merged_frame.stack())
        merged_frame.reset_index(inplace=True)
        merged_frame.columns = (
            ["ext_gene", "ens_gene", "target_id", "Sample", col]
            if sleuth else
            ["target_id", "Sample", col]
        )
        merged_frame.set_index(
            (["ext_gene", "ens_gene", "target_id"] if sleuth else "target_id"),
            inplace=True)

    if (tr2g is not None) and (not sleuth):
        logging.debug("Adding Genes names")
        merged_frame = pd.merge(
            merged_frame,
            read_t2g(tr2g),
            left_index=True,
            right_index=True,
            how="left"
        )

    if stack or drop_null or drop_na:
        logging.debug("Merged frame post processed.")
        logging.debug(merged_frame.shape)
        logging.debug(merged_frame.head())

    logging.debug("Writing results")
    merged_frame.to_csv(
        str(destin),
        sep="\t",
        index=True,
        header=True,
        index_label=("target_id" if index_label else False)
    )
