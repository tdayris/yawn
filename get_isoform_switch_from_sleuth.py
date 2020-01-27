#!/usr/bin/python3
# conding: utf-8

import begin
import logging
import pandas as pd
import traceback


def filtr_switch(files):
    yield from [
        (
            f[f["ens_gene"].isin(f[f[stat] <= threashold]["ens_gene"])],
            p,
            stat,
            threashold
        )
        for f, p in [
            (pd.read_csv(f, sep="\t", index_col=0, header=0), f)
            for f in files
        ]
        for stat in ["qval", "pval"]
        for threashold in [0.001, 0.005, 0.01, 0.05]
    ]


@begin.start
@begin.logging
@begin.tracebacks
@begin.convert(diffexp=begin.utils.tobool)
def main(*files: "Path to multiple Sleuth files"):
    for r, p, s, t in (filtr_de(files) if diffexp else filtr_switch(files)):
        logging.debug("Working on %s, with %s: %s" % (p, s, str(t)))
        r.to_csv("%s_isoform_switch_%s_%s.tsv" % (p, s, str(t)), sep="\t")
