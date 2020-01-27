#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");
base::library("BSgenome.Hsapiens.UCSC.hg38");

orf_analysed <- base::readRDS(
  base::as.character(snakemake@input[["orf"]])
);

orf_analysed <- IsoformSwitchAnalyzeR::extractSequence(
  orf_analysed,
  genomeObject = BSgenome.Hsapiens.UCSC.hg38,
  removeORFwithStop = FALSE,
  writeToFile = TRUE,
  pathToOutput = file.path(getwd(), "ISAR"),
  outputPrefix = snakemake@params[["prefix"]]
);

base::saveRDS(
  orf_analysed,
  base::as.character(snakemake@output[["sequence"]])
);
