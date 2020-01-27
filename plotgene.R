#!/usr/bin/R

library("IsoformSwitchAnalyzeR");

message("Libraries loaded");

subset_switch <- base::readRDS(
  base::as.character(snakemake@input[["switches"]])
);

nb_switches <- dim(
  subset_switch$isoformFeatures[subset_switch$isoformFeatures$isoform_switch_q_value < 0.05, ]
)[1];

message("Switches loaded");

switches <- IsoformSwitchAnalyzeR::extractTopSwitches(
  subset_switch,          # The switch analysis
  n = nb_switches,        # Number of switching features
  sortByQvals = TRUE,     # Sort by qvals
  extractGenes = FALSE    # Estract all isoforms (not genes only)
);

message("Switches extracted");

write.table(switches, file = snakemake@output[["results"]]);
message("Table saved");

for (i in seq(1, dim(switches)[1])) {
  print(switches[i, ]);
  gene_name <- base::as.character(switches[i,]$gene_name);
  gene_id <- base::as.character(switches[i,]$gene_id);
  cond1 <- base::as.character(switches[i,]$condition_1);
  cond2 <- base::as.character(switches[i,]$condition_2);
  name <- base::paste0(
    gene_name, "_", gene_id,
    "_",
    cond1, "_vs_", cond2, ".png"
  );
  png_path <- base::file.path(snakemake@output[["images"]], name);
  message(paste("Plot", name));
  png(png_path, units = "px", height = 934, width = 1245);
  IsoformSwitchAnalyzeR::switchPlot(
    subset_switch, gene = gene_id, condition1 = cond1, condition2 = cond2
  );
  dev.off();
}
message("All png saved");
