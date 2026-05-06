#!/usr/bin/env Rscript
# Build Reijns rosetta CSVs and ID476 catalogs.
# Output files are written next to this script in inst/.

library(mSigSpectra)

ca <- commandArgs(trailingOnly = FALSE)
m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
inst_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()

source(file.path(inst_dir, "make_reijns_rosetta_csv.R"))
source(file.path(inst_dir, "make_reijns_catalog.R"))

make_reijns_rosetta_csv(
  vcf_dir    = inst_dir,
  ref_genome = "hg38",
  out_path   = file.path(inst_dir, "reijns_cells_keep_high_R.csv"),
  drop_high_R = FALSE
)

make_reijns_rosetta_csv(
  vcf_dir    = inst_dir,
  ref_genome = "hg38",
  out_path   = file.path(inst_dir, "reijns_cells_drop_high_R.csv"),
  drop_high_R = TRUE
)

make_reijns_catalog(
  vcf_dir  = inst_dir,
  ref_genome = "hg38",
  out_keep = file.path(inst_dir, "reijns_cells_catalog_keep_high_R.tsv"),
  out_drop = file.path(inst_dir, "reijns_cells_catalog_drop_high_R.tsv")
)
