#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Intention:
- Compute fold change from a count tsv-formatted file.

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

import argparse
import begin
import itertools
import logging
import numpy as np
import os
import pandas as pd
import traceback


@begin.start(formatter_class=argparse.RawTextHelpFormatter)
@begin.logging
@begin.tracebacks
@begin.convert(numreads=str, condition_array=str, output=str,
               output_excel=begin.utils.tobool)
def main(numreads: "Path to a tsv file with the raw number of reads",
         *condition_array: "Space separated list of conditions",
         output: "Path to output file"="FC.tsv"):
    """
    Read a NumRead file from a tsv-formatted count-file and return a FoldChange
    This script accepts n conditions and will perform permutations over the
    given conditions in order to produce all FC.

    The FoldChange is defined as :
    * -1 / ratio if ration < 1
    * ratio if ratio is => 1
    * 'inf' if denominator is 0

    ratio = mean(counts obs) / mean(count ref)
    """
    foldchanges = pd.DataFrame()
    # Loading datasets
    if os.path.exists(output):
        raise FileExistsError("%s already exists" % output)

    data = pd.read_csv(numreads, header=0, sep="\t", index_col=0)
    data = data[list(data.select_dtypes(include=[np.number]).columns.values)]
    samples = data.columns.tolist()
    logging.debug("Data loaded")

    if not len(condition_array) == len(samples):
        raise ValueError("Uncoherent number of conditions and samples")

    # Computing ratio
    logging.debug("Processing Fold Change")
    samples_to_conditions = {}
    for s, c in zip(samples, condition_array):
        try:
            samples_to_conditions[c].append(s)
        except KeyError:
            samples_to_conditions[c] = [s]

    logging.debug("The mapping sample/condition results "
                  "in: %s" % str(samples_to_conditions))

    for c, s in samples_to_conditions.items():
        logging.debug("Computing mean count per target over %s" % c)
        data["mean(%s)" % c] = data[s].mean(axis=1)

    for cs in itertools.permutations(samples_to_conditions.keys(), 2):
        logging.debug("Computing FC for %s / %s" % cs)
        c1, c2 = cs
        ratio = data["mean(%s)" % c1] / data["mean(%s)" % c2]
        # ratio = data["mean(%s)" % c1].div(data["mean(%s)" % c2])
        foldchanges["FC(%s/%s)" % cs] = [
            -1/r if 0 < r < 1 else (pd.np.NaN if r == 0 else r) for r in ratio
        ]
    foldchanges["GeneIdentifier"] = data.index.tolist()
    foldchanges.set_index("GeneIdentifier", inplace=True)

    # Saving
    logging.debug("Saving data")
    foldchanges.to_csv(output, sep="\t")
