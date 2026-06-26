#!/usr/bin/env Rscript
# Plot SBS96 and ID83 catalogs in this folder using mSigPlot.

suppressPackageStartupMessages({
  pkgload::load_all("/home/steve/github/mSigSpectra")
  library(mSigPlot)
})

this_dir <- "/home/steve/github/mSigSpectra/development/duocarmycin-cell-line-exposures"

read_and_bind <- function(files) {
  cats <- lapply(files, function(f) {
    cat <- read_catalog(f, region = "genome")
    colnames(cat) <- sub("\\.(SBS96|ID83|DBS78)\\.csv$", "", basename(f))
    cat
  })
  cbind_catalogs(cats)
}

read_simple_matrix <- function(files, suffix) {
  cols <- lapply(files, function(f) {
    df <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    v <- df[[2L]]; names(v) <- df$MutationType
    v
  })
  m <- do.call(cbind, cols)
  colnames(m) <- sub(paste0("\\.", suffix, "\\.csv$"), "", basename(files))
  m
}

sbs_files   <- sort(list.files(this_dir, pattern = "\\.SBS96\\.csv$", full.names = TRUE))
id83_files  <- sort(list.files(this_dir, pattern = "\\.ID83\\.csv$",  full.names = TRUE))
id89_files  <- sort(list.files(this_dir, pattern = "\\.ID89\\.csv$",  full.names = TRUE))
id476_files <- sort(list.files(this_dir, pattern = "\\.ID476\\.csv$", full.names = TRUE))

plot_SBS96_pdf(read_and_bind(sbs_files),
               filename = file.path(this_dir, "duocarmycin_SBS96.pdf"))
plot_ID83_pdf(read_and_bind(id83_files),
              filename = file.path(this_dir, "duocarmycin_ID83.pdf"))
plot_ID89_pdf(read_simple_matrix(id89_files, "ID89"),
              filename = file.path(this_dir, "duocarmycin_ID89.pdf"))
plot_ID476_pdf(read_simple_matrix(id476_files, "ID476"),
               filename = file.path(this_dir, "duocarmycin_ID476.pdf"))

message("Wrote SBS96, ID83, ID89, ID476 PDFs in ", this_dir)
