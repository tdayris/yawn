#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");

splice_analysed <- base::readRDS(
  base::as.character(snakemake@input[["splice"]])
);

orf_analysed <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
  switchAnalyzeRlist = splice_analysed,
);

base::saveRDS(
  orf_analysed,
  base::as.character(snakemake@output[["consequences"]])
);
