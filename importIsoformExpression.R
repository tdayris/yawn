#!/usr/bin/R

base::library("IsoformSwitchAnalyzeR");

isar_input <- IsoformSwitchAnalyzeR::importIsoformExpression(
  parentDir = base::as.character(snakemake@params[["parent_dir"]]),
  addIsofomIdAsColumn = TRUE
);

message("Expression loaded");

design <- utils::read.table(
  base::as.character(snakemake@params[["design"]]),
  sep = "\t",
  header = TRUE
);

design <- design[c("Sample_id", snakemake@params[["condition"]])];
colnames(design) <- c("sampleID", "condition");
print(head(design));
message("Table red.");

switch_list <- IsoformSwitchAnalyzeR::importRdata(
  isoformCountMatrix = isar_input$counts,
  isoformRepExpression = isar_input$abundance,
  designMatrix = design,
  isoformExonAnnoation = as.character(snakemake@params[["gtf"]]),
  isoformNtFasta = as.character(snakemake@params[["transcripts"]]),
  addAnnotatedORFs = TRUE
);

message("Switch list created");

filtered <- IsoformSwitchAnalyzeR::preFilter(
  switch_list,
  geneExpressionCutoff = 0.1,
  removeSingleIsoformGenes = TRUE,
);

message("Filters applied");

base::saveRDS(
  filtered,
  base::as.character(snakemake@output[["rds"]])
);
