#!/usr/bin/Rscript

# Loading libraries
base::suppressMessages({
  base::library(maftools);  # Handle VCF and MAF
  base::library(argparse);  # Parse command line
});
base::message("Libraries loaded");

# Command line parsing
parser <- argparse::ArgumentParser(
  description = "Perform MAFTools graphics on your Annovar VCF file"
);

parser$add_argument(
  "VCF",
  metavar = "PATH",
  help    = "Path to your VCF file",
  type    = "character",
  nargs   = "+"
);

parser$add_argument(
  "MAF",
  metavar = "PATH",
  help    = "Prefix to output file",
  type    = "character"
);

parser$add_argument(
  "--genes",
  metavar = "GENE",
  help    = "Space separated list of genes of interest",
  type    = "character",
  nargs   = "+",
  default = c("TP53")
)

parser$add_argument(
  "--ncbi_build",
  help    = "NCBI_Build field in MAF file will be filled with this value.",
  type    = "character",
  default = "hg19"
);

parser$add_argument(
  "--table",
  help    = "reference table used for gene-based annotations.",
  type    = "character",
  default = "ensGene"
);

opt <- parser$parse_args();
base::message("Argument parsed:");
print(opt)

message("Working on the following VCF files:");
print(opt$VCF);

# File conversion
maf_file <- maftools::annovarToMaf(
  annovar  = opt$VCF,
  refBuild = opt$ncbi_build,
  table    = opt$table,
  ens2hugo = TRUE,
  basename = opt$MAF,
  MAFobj   = TRUE
);



# Plotting data
grDevices::png(
  filename = base::paste0(opt$MAF, "_Summary.png"),
  width    = 1024,
  height   = 1024
);
maftools::plotmafSummary(
  maf       = maf_file,
  rmOutlier = TRUE,
  addStat   = 'median',
  dashboard = TRUE,
  titvRaw   = FALSE
);
grDevices::dev.off();

plot_gene <- function(gene) {
  message(
    paste(
      "Plotting",
      gene,
      sep = ": "
    )
  );
  grDevices::png(
    filename = base::paste0(opt$MAF, "_", gene, ".png"),
    width    = 1024,
    height   = 1024
  );
  maftools::lollipopPlot(
    maf             = maf_file,
    gene            = gene,
    AACol           = 'Gene.refGene',
    showDomainLabel = TRUE
  );
  grDevices::dev.off();
}

sapply(opt$genes, function(gene_name) plot_gene(gene_name))

grDevices::png(
  filename = base::paste0(opt$MAF, "_Top_10_Oncoplot.png"),
  width    = 1024,
  height  = 1024
);
maftools::oncoplot(
  maf = maf_file,
  top = 10
);
grDevices::dev.off();

grDevices::png(
  filename = base::paste0(opt$MAF, "_Interest_Oncoplot.png"),
  width    = 1024,
  height   = 1024
);
maftools::oncoplot(
  maf   = maf_file,
  genes = opt$genes
);
grDevices::dev.off();

maf_titv <- maftools::titv(
  maf    = maf_file,
  plot   = FALSE,
  useSyn = TRUE
);
grDevices::png(
  filename = base::paste0(opt$MAF, "_TITV.png"),
  width    = 1024,
  height   = 1024
);
maftools::plotTiTv(
  res = maf_titv
);
grDevices::dev.off();
