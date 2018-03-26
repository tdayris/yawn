#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
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


def read_salmons(prefix, suffix, *files) -> pd.DataFrame:
    """
    Returns a unique DataFrame from multiple Salmon files. It contains:
        - TPM for each sample
        - Mean TPM
    """
    data = None
    for salmon in files:
        logging.debug("Parsing: %s" % salmon)
        tmp = pd.read_csv(salmon, sep="\t", header=0, index_col=0)[["TPM"]]
        sample_id = str(salmon)
        sample_id = sample_id if suffix == "" else sample_id[:-len(suffix)]
        sample_id = sample_id if prefix == "" else sample_id[len(prefix):]
        tmp.columns = [sample_id]
        logging.debug(tmp.head())
        try:
            data = pd.merge(data, tmp, left_index=True, right_index=True)
        except ValueError:
            data = tmp

    data["MeanTPM"] = data.mean(axis=1)
    return data


@begin.start
@begin.logging
@begin.tracebacks
def main(condition1: "Path to files from condition 1 (comma separated)",
         condition2: "Path to files from condition 2 (comma separated)",
         condition1_name: "The condition 1 name"="Condition1",
         condition2_name: "The condition 2 name"="Condition2",
         output: "Path to output file"="FoldChange.tsv",
         prefix: "The prefix to be removed to get "
                 "sample name from its file name"="",
         suffix: "The suffix to be removed to get "
                 "sample name from its file name"="",) -> None:
    """
    This script reads salmon files and returns a table of foldchanges
    """
    print(condition1.split(","))
    cond1 = read_salmons(prefix, suffix, *condition1.split(","))
    cond2 = read_salmons(prefix, suffix, *condition2.split(","))
    suf1 = "_%s" % condition1_name
    suf2 = "_%s" % condition2_name
    logging.debug("Mergind data")
    fc = pd.merge(
        cond1, cond2,
        left_index=True, right_index=True,
        suffixes=[suf1, suf2],
        how="outer"
    )
    fc.fillna(0)
    logging.debug("Conputing Fold Change")
    fc_name = "FoldChange %s/%s" % (condition1_name, condition2_name)
    fc[fc_name] = fc["MeanTPM%s" % suf1] / fc["MeanTPM%s" % suf2]
    logging.debug(fc.head())
    logging.debug("Saving to %s" % output)
    fc.to_csv(output, sep="\t")
