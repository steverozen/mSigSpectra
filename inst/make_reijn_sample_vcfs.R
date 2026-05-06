#!/usr/bin/env Rscript
# make_reijn_sample_vcfs.R
#
# Split Reijn.rpe1.vcfs.with.annotation.txt into one per-sample VCF file,
# keeping only the standard VCF columns (CHROM through FORMAT).
# Output files are named reijn_mut_cell_<Sample>.vcf.

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(data.table))

ca       <- commandArgs(trailingOnly = FALSE)
m        <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
inst_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()

p <- arg_parser("Split Reijn RPE1 annotation file into per-sample VCFs")
p <- add_argument(p, "--input",
  default = "/home/steve/github/Code_Liu_2025/Manuscript_data/Reijn.rpe1.vcfs.with.annotation.txt",
  help    = "Path to tab-separated Reijn annotation file")
p <- add_argument(p, "--out-dir",
  default = inst_dir,
  help    = "Directory for output .vcf files")
args <- parse_args(p)

message("Reading ", basename(args$input))
dt <- data.table::fread(args$input, sep = "\t", header = TRUE, quote = "")
message("  ", nrow(dt), " rows, ", ncol(dt), " columns")

vcf_cols <- colnames(dt)[1:9]  # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT

for (s in unique(dt$Sample)) {
  sub <- dt[Sample == s, ..vcf_cols]
  out <- file.path(args$out_dir, paste0("reijn_cell_", s, ".vcf"))
  writeLines(paste0("#", paste(vcf_cols, collapse = "\t")), out)
  data.table::fwrite(sub, out, sep = "\t", col.names = FALSE,
                     append = TRUE, quote = FALSE)
  message("Wrote ", out, " (", nrow(sub), " rows)")
}
