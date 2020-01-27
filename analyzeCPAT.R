#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");

sequence_analysed <- base::readRDS(
  base::as.character(snakemake@input[["sequence"]])
);

orf_analysed <- IsoformSwitchAnalyzeR::analyzeCPAT(
  switchAnalyzeRlist = sequence_analysed,
  pathToCPATresultFile = base::as.character(snakemake@input[["cpat"]]),
  codingCutoff = 0.725,
  removeNoncodinORFs = TRUE
);

base::saveRDS(
  orf_analysed,
  base::as.character(snakemake@output[["cpat"]])
);
