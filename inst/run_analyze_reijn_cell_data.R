#!/usr/bin/env Rscript
# Run make_reijn_rosetta_csv() twice: once keeping R > 9 rows, once dropping them.
# Output CSVs are written next to this script in inst/.

library(mSigSpectra)

ca <- commandArgs(trailingOnly = FALSE)
m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
inst_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()

source(file.path(inst_dir, "make_reijn_rosetta_csv.R"))

make_reijn_rosetta_csv(
  vcf_dir = inst_dir,
  ref_genome = "hg38",
  out_path = file.path(inst_dir, "reijn_cells_keep_high_R.csv"),
  drop_high_R = FALSE
)

make_reijn_rosetta_csv(
  vcf_dir = inst_dir,
  ref_genome = "hg38",
  out_path = file.path(inst_dir, "reijn_cells_drop_high_R.csv"),
  drop_high_R = TRUE
)
