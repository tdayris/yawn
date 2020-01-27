#!/usr/bin/python3

import begin
import logging
import pandas as pd
import traceback


@begin.start
@begin.convert(tsv=str, tr2g=str, transcripts=begin.utils.tobool)
@begin.logging
@begin.tracebacks
def main(*tsv: "Path to multiple TSV files with genes ID in first column",
         tr2g: "Path to a transcript to gene table",
         transcripts: "The id in first column are transcript id"
                      " and not gene ones" = False):
    """
    This script adds a final column with genes names (symbol from HUGO
    database) in a TSV file. The genes ID must be in first column!
    """
    print(tr2g)
    t2g = pd.read_csv(str(tr2g), header=None, index_col=1, sep="\t")
    t2g.columns = ["Ensembl_ID", "Hugo_ID"]
    logging.debug(t2g.head())

    if not transcripts:
        t2g = (t2g.reset_index()[["Ensembl_ID", "Hugo_ID"]]
                  .drop_duplicates()
                  .set_index("Ensembl_ID"))
        logging.debug(t2g.head())

    for tsv_file in tsv:
        logging.debug("Working on %s" % tsv_file)
        data = pd.merge(
            pd.read_csv(tsv_file, header=0, index_col=0, sep="\t"),
            t2g,
            left_index=True,
            right_index=True,
            how="left"
        )
        logging.debug(data.head())

        logging.debug("Writing results for: %s" % tsv_file)
        logging.debug("Writing to: %s" % "%sannotated.tsv" % tsv_file[:-3])
        data.to_csv(
            "%sannotated.tsv" % tsv_file[:-3],
            sep="\t",
            index=True,
            header=True
        )
