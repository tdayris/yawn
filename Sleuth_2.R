#!/usr/bin/Rscript
# Author: Dayris Thibault

# Intention:
# - Perform a differential analysis with Sleuth
#
# Unlincence terms of use:
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>
#
# Edition guidelines:
# version: [core].[version].[edition]
# namespace: use_the_snake_case for var names, UseCamelCase for class names
# rmemeber: _function is private by intention
#
# Add your name in author list if you modify anything. Cheers and thank you!


suppressMessages({
  library("sleuth");
  library("argparse");
  library("dplyr");
  library("wasabi");
});
message("Libraries loaded");

parser <- argparse::ArgumentParser(
  description = "Perform your Sleuth differential analysis"
);

parser$add_argument(
  "design",
  metavar = "PATH",
  help = "Path to the design file",
  type = "character"
);

parser$add_argument(
  "formula",
  metavar = "FORMULA",
  help = "The R formula to be used to perform the analysis",
  type = "character"
);

parser$add_argument(
  "intercept",
  metavar = "FORMULA",
  help = "The intercept model",
  type = "character"
);

parser$add_argument(
  "-w", "--workdir",
  metavar = "PATH",
  default = getwd(),
  help = "Path to the working directory",
  type = "character"
);

parser$add_argument(
  "-t", "--tr2gene",
  metavar = "PATH",
  default = NULL,
  help = "Path to the transcript to gene table",
  type = "character"
);

parser$add_argument(
  "--qval",
  metavar = "FLOAT",
  default = 0.05,
  help = "Q-Value threshold to filter out genes",
  type = "double"
)

parser$add_argument(
  "-p", "--process",
  metavar = "INT",
  default = 1,
  help = "Maximum number of threads",
  type = "integer"
);

parser$add_argument(
  "-o", "--output",
  metavar = "NAME",
  default = "Sleuth_Results",
  help = "Base name of your sleuth results",
  type = "character"
);

parser$add_argument(
  "-a", "--genes_analysis",
  # action = "store_true",
  help = "Perform a gene analysis instead of a transcript one",
  default = FALSE,
  type = "logical"
);
message("ArgumentParser defined");

# Parse arguments
opt <- parser$parse_args();
message("Argument parsed");

opt$design <- as.character(opt$design);
opt$formula <- as.formula(opt$formula);
opt$intercept <- as.formula(opt$intercept);
message("Types re-defined");

threads <- 1
mc.cores <- threads;
options(mc.cores = threads);
message("Multi threading options set");


# Define output files names
opt$output_table <- file.path(
  opt$workdir,
  paste0(opt$output)
);
opt$output_rds <- file.path(
  opt$workdir,
  paste0(opt$output, ".rds")
);
opt$output_exp <- file.path(
  opt$workdir,
  paste0(opt$output, ".expression.tsv")
);
opt$args_out <- file.path(
  opt$workdir,
  paste0(opt$output, ".arguments.rds")
);
message("Output paths defined");

# Read transcripts to gene table if any
if (! is.null(opt$tr2gene)) {
  opt$tr2gene <- as.character(opt$tr2gene);
  opt$t2g <- read.table(
    opt$tr2gene,
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t"
  );
  colnames(opt$t2g) <- c("ens_gene", "target_id", "ext_gene");
} else {
  opt$t2g <- NULL
}
message("Transcript to Gene table defined (if needed)");

# Set tmp directory to avoid /tmp issues
TMP <- file.path(opt$workdir, "tmp");
TEMP <- file.path(opt$workdir, "tmp");
TMPDIR <- file.path(opt$workdir, "tmp");
message("TMP directories set");

# Load meta-data table
opt$meta <- read.table(opt$design, header = TRUE, sep = "\t");
colnames(opt$meta)[1:2] <- c("sample", "path");
opt$meta$path <- as.character(opt$meta$path);
message("Meta data loaded");

# Convert from salmon to Kallisto
sapply(opt$meta$path, function(path) wasabi::prepare_fish_for_sleuth(path));
message("Converted with Wasabi");

# Save arguments to help quick re-run and complements
saveRDS(opt, file = opt$args_out);
message("Arguments parsed");

# Preparing sleuth analysis
if (opt$genes_analysis) {
    so <- sleuth::sleuth_prep(
        sample_to_covariates = opt$meta,
        full_model = as.formula(opt$formula),
        target_mapping = opt$t2g,
        aggregation_column = "ens_gene",
        extra_bootstrap_summary = TRUE,
        num_cores = threads
    );
} else {
    so <- sleuth::sleuth_prep(
        sample_to_covariates = opt$meta,
        full_model = as.formula(opt$formula),
        target_mapping = opt$t2g,
        extra_bootstrap_summary = TRUE,
        num_cores = threads
    );
}
message("Sleuth object prepared");

# Starting analysis
so <- sleuth::sleuth_fit(
  obj = so,
  formula = as.formula(opt$formula),
  fit_name = "full"
);

so <- sleuth::sleuth_fit(
  obj = so,
  formula = as.formula(opt$intercept),
  fit_name = "reduced"
);
message("Sleuth model fitted");

so <- sleuth::sleuth_lrt(so, "reduced", "full");
message("Performed likehood ratio test");

for (beta in colnames(so$design_matrix)) {
  message(paste0("Performed Wald Test on ", beta));
  so <- sleuth_wt(obj = so, which_model = "full", which_beta = beta);
}
message("Wald tests over.");

saveRDS(so, opt$output_rds);
message("RDS object saved.");

# Save multiple results per beta
for (beta in colnames(so$design_matrix)){
  # Get result table
  betatable <- sleuth::sleuth_results(
    obj = so,
    test = beta,
    which_model = "full",
    test_type = "wt"
  );

  # Save it!
  write.table(
    betatable,
    file = paste0(opt$output_table, "_", beta, ".tsv"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  );
  message(paste("Beta table", beta, "saved."));

  # Extract only the columns that will please the biologists
  # if (is.null(opt$gene_analysis)) {
  #   betatable <- betatable %>%
  #     dplyr::select(ext_gene, target_id, qval, b);
  # } else {
  #   betatable <- betatable %>%
  #     dplyr::select(ext_gene, target_id, qval, b, ens_gene);
  # }
  # write.table(
  #   betatable,
  #   file = paste0(opt$output_table, "_", beta, "_readable.tsv"),
  #   quote = FALSE,
  #   row.names = FALSE,
  #   sep = "\t"
  # );
  # message(paste("Readable beta table", beta, "saved."));

  # Rename columns and filter the table for GSEApp
  if (is.null(opt$gene_analysis)) {
    betatable <- dplyr::filter(betatable, qval <= opt$qval) %>%
      setNames(c(
        "GeneIdentifier", "EnsemblTranscriptID", "Q-Value", "stat_change"
    ));
  } else {
    betatable <- dplyr::filter(betatable, qval <= opt$qval) %>%
      setNames(c(
        "GeneIdentifier", "EnsemblTranscriptID", "Q-Value",
        "stat_change", "EnsemblGeneIdentifier"
    ));
  }
  write.table(
    betatable,
    file = paste0(opt$output_table, "_", beta, "_gseapp.tsv"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  );
  message(paste("GSEApp compatible beta table", beta, "saved."));
}
message("Sleuth Result table saved");

# Geathering transcripts expression
quant <- sleuth::kallisto_table(so) %>%
  dplyr::select(target_id, sample, tpm) %>%
  tidyr::spread(sample, tpm);

write.table(
  quant,
  file = opt$output_exp,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
);
message("Transcripts expression table saved");

message("Everything is done!");
