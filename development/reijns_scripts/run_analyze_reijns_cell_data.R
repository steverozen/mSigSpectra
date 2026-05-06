#!/usr/bin/env Rscript
# Build Reijns rosetta CSVs and ID476 catalogs.
# Reads VCFs from <repo>/inst/reijns_cell_data/.
# Output CSV / TSV files are written next to this script.

library(mSigSpectra)

ca <- commandArgs(trailingOnly = FALSE)
m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
script_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()
vcf_dir <- normalizePath(
  file.path(script_dir, "..", "..", "inst", "reijns_cell_data")
)

source(file.path(script_dir, "make_reijns_rosetta_csv.R"))
source(file.path(script_dir, "make_reijns_catalog.R"))

make_reijns_rosetta_csv(
  vcf_dir    = vcf_dir,
  ref_genome = "hg38",
  out_path   = file.path(script_dir, "reijns_cells_keep_high_R.csv"),
  drop_high_R = FALSE
)

make_reijns_rosetta_csv(
  vcf_dir    = vcf_dir,
  ref_genome = "hg38",
  out_path   = file.path(script_dir, "reijns_cells_drop_high_R.csv"),
  drop_high_R = TRUE
)

make_reijns_catalog(
  vcf_dir  = vcf_dir,
  ref_genome = "hg38",
  out_keep = file.path(script_dir, "reijns_cells_catalog_keep_high_R.tsv"),
  out_drop = file.path(script_dir, "reijns_cells_catalog_drop_high_R.tsv")
)
