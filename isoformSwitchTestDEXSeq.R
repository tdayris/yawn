#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");

filtered <- base::readRDS(
  base::as.character(snakemake@input[["rds"]])
);

switches_analysed <- IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = filtered,
  reduceToSwitchingGenes = TRUE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  alpha = 0.05
);

extractSwitchSummary(
  switches_analysed
);

base::saveRDS(
  switches_analysed,
  base::as.character(snakemake@output[["switches"]])
);
