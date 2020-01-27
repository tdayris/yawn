#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");

pfam_analyzed <- base::readRDS(
  base::as.character(snakemake@input[["pfam"]])
);

splice_analysed <- IsoformSwitchAnalyzeR::analyzeAlternativeSplicing(
  switchAnalyzeRlist = pfam_analyzed,
  onlySwitchingGenes = FALSE,
);

base::saveRDS(
  splice_analysed,
  base::as.character(snakemake@output[["splice"]])
);
