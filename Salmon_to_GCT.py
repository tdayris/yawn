#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script builds a GCT file for GSEA, from multiple Salmon files.

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
import os.path
import pandas as pd
import traceback

from typing import List


def read_Salmon(prefix=str, suffix=str,
                *salmon_files: List[str]) -> pd.DataFrame:
    """
    Read multiple Salmon files and return a DataFrame
    """
    merged_results = None
    for salmon in salmon_files:
        logging.debug("Working on %s" % salmon)
        tmp = pd.read_csv(
            salmon,
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
        tmp = tmp[["TPM"]]
        sample_id = salmon
        sample_id = sample_id if suffix == "" else sample_id[:-len(suffix)]
        sample_id = sample_id if prefix == "" else sample_id[len(prefix):]
        logging.debug("%s's shape: %s" % (sample_id, str(tmp.shape)))
        tmp.columns = [sample_id]
        try:
            merged_results = pd.merge(
                merged_results,
                tmp,
                left_index=True,
                right_index=True
            )
        except ValueError:
            merged_results = tmp

    logging.debug("Frames merged:")
    return merged_results


@begin.start
@begin.convert(output_file=str)
@begin.logging
@begin.tracebacks
def main(*salmon_files: "Space separated list of salmon files List[str]",
         t2g_file: "Path to a T2G file",
         prefix: "The prefix to be removed to get "
                 "sample name from its file name"="",
         suffix: "The suffix to be removed to get "
                 "sample name from its file name"="",
         metadata_file: "Path to a metadata file"=None,
         output_file: "Path to output file"="salmon.gct") -> None:
    """
    This script builds a GCT file for GSEA, from multiple Salmon files.
    """
    if os.path.exists(output_file):
        raise FileExistsError("%s already exists" % output_file)

    salmon_frame = read_Salmon(prefix, suffix, *salmon_files)
    # Add the required salmon description column
    if metadata_file is None:
        salmon_frame["Description"] = [None]*len(salmon_frame)
    else:
        logging.debug("Reading %s" % metadata_file)
        metadata = pd.read_csv(
            metadata_file, sep="\t", index_col=0, header=None
        )

        salmon_frame = pd.merge(
            salmon_frame,
            metadata,
            left_index=True,
            right_index=True
        )

        salmon_frame.columns = (
            salmon_frame.columns.tolist()[0:-1] + ["Description"]
        )

    # Replacing ensembl ID with HUGO ID
    logging.debug("Replacing Ensembl ID by Hugo names")
    t2g = pd.read_csv(t2g_file, sep="\t", header=None, index_col=1)
    t2g.columns = ["Gene_ID", "Name"]
    t2g = t2g.apply(lambda x: x.astype(str).str.upper())

    salmon_frame = pd.merge(
        salmon_frame,
        t2g,
        left_index=True,
        right_index=True
    )

    salmon_frame.reset_index(inplace=True)
    salmon_frame.set_index("Name", inplace=True)
    del salmon_frame["Gene_ID"]
    del salmon_frame["index"]

    # Reorder the columns as expected
    salmon_cols = salmon_frame.columns.tolist()
    salmon_frame = salmon_frame[
        ["Description", *salmon_cols[0:-1]]
    ]

    # Save file with comment lines
    nline, ncol = salmon_frame.shape
    ncol -= 1
    logging.debug("Saving file: %s" % output_file)
    with open(output_file, "w") as outfile:
        outfile.write("#1.2\n%i\t%i\n" % (nline, ncol))
        salmon_frame.index.str.upper()
        salmon_frame.to_csv(outfile, sep="\t")
