#!/usr/bin/env Rscript
# Run make_rosetta_csv() twice: once dropping rows with R > 9, once keeping them.
# Output CSVs are written next to this script in the same directory.

library(mSigSpectra)

ca       <- commandArgs(trailingOnly = FALSE)
m        <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
inst_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()

source(file.path(inst_dir, "make_rosetta_csv.R"))

make_rosetta_csv(
  vcf_dir     = inst_dir,
  ref_genome  = "GRCh37",
  out_path    = file.path(inst_dir, "rosetta_drop_high_R.csv"),
  drop_high_R = TRUE
)

make_rosetta_csv(
  vcf_dir     = inst_dir,
  ref_genome  = "GRCh37",
  out_path    = file.path(inst_dir, "rosetta_keep_high_R.csv"),
  drop_high_R = FALSE
)
