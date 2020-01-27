#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");

sequence_analysed <- base::readRDS(
  base::as.character(snakemake@input[["cpat"]])
);

orf_analysed <- IsoformSwitchAnalyzeR::analyzePFAM(
  switchAnalyzeRlist = sequence_analysed,
  pathToPFAMresultFile = base::as.character(snakemake@input[["pfam"]])
);

base::saveRDS(
  orf_analysed,
  base::as.character(snakemake@output[["pfam"]])
);
