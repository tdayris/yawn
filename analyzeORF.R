#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");
base::library("BSgenome.Hsapiens.UCSC.hg38");

switches_analysed <- base::readRDS(
  base::as.character(snakemake@input[["switches"]])
);

orf_analysed <- IsoformSwitchAnalyzeR::analyzeORF(
  switchAnalyzeRlist = switches_analysed,
  genomeObject = BSgenome.Hsapiens.UCSC.hg38,
  orfMethod = "longest"
);

base::saveRDS(
  orf_analysed,
  base::as.character(snakemake@output[["orf"]])
);
