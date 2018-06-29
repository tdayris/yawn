#!/usr/bin/python3
# conding: utf-8

import begin
import logging
import pandas as pd
import traceback


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(sleuth_files=str, p_threashold=float, q_threashold=float)
def main(*sleuth_files: "Path to multiple Sleuth files",
         p_threashold: "P-Value threashold"=0.001,
         q_threashold: "Q-Value threashold"=0.001,):
    """
    This script subsets the Sleuth files
    """
    for file in sleuth_files:
        data = pd.read_csv(file, sep="\t", index_col=0, header=0)
        data = data[
            (data["qval"] <= q_threashold) & (data["pval"] <= p_threashold)
        ]
        data = data[["ens_gene", "target_id", "pval", "qval", "b"]]
        if len(data) >= 0:
            data.to_csv(
                "%s_qval%s_pval%s.tsv" % (file, q_threashold, p_threashold),
                sep="\t"
            )
        else:
            logging.warning(
                "No data available for %s, with %f" % (file, q_threashold)
            )
