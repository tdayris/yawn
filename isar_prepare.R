#!/usr/bin/R

# Prepare the ISA-R object based on input counts and
# sleuth idfferential analysis.


# Loading libraries
suppressMessages({
  base::library("sleuth");
  base::library("argparse");
  base::library("dplyr");
  base::library("IsoformSwitchAnalyzeR");
});
message("Libraries loaded");

# Parsing command line
parser <- argparse::ArgumentParser(
  description = "Prepare your Isoform Swith Analysis with Sleuth results and Salmon counts"
);

parser$add_argument(
  "salmon_quants",
  metavar = "salmon",
  help = "Path to the salmon quantification directory",
  type = "character"
);

parser$add_argument(
  "-t", "--threads",
  help = "Maximum number of threads used",
  type = "integer",
  default = 1
);


message("ArgumentParser defined");

# Parse arguments
opt <- parser$parse_args();
message("Argument parsed");

# Casting variables
opt$salmon_quants <- as.character(opt$salmon_quants);
message("Types re-defined");

# Threading behaviour
threads <- 1
mc.cores <- threads;
options(mc.cores = threads);
message("Multi threading options set");

# Set tmp directory to avoid /tmp issues
TMP <- file.path(opt$workdir, "tmp");
TEMP <- file.path(opt$workdir, "tmp");
TMPDIR <- file.path(opt$workdir, "tmp");
message("TMP directories set");

