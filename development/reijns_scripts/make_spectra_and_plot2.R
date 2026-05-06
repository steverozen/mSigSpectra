# make_spectra_and_plot2.R — one-off
#
# Read every indel VCF in this folder, build ID83, ID89, and ID476
# spectra catalogs with mSigSpectra, and write one multi-page PDF of
# mSigPlot plots next to the VCFs. Each page holds one sample's three
# panels using the fig1.R layout (see ../one_page_83_89_476.R).

library(mSigSpectra)
library(mSigPlot)
library(ggplot2)

# Locate the directory that contains this script so siblings resolve whether
# sourced interactively, via Rscript, or via source(chdir=).
this_script <- tryCatch(
  {
    ca <- commandArgs(trailingOnly = FALSE)
    m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
    if (length(m)) {
      m[1]
    } else {
      ofile <- NULL
      for (i in seq_len(sys.nframe())) {
        f <- sys.frame(i)
        if (!is.null(f$ofile) && is.character(f$ofile) && nzchar(f$ofile)) {
          ofile <- f$ofile
          break
        }
      }
      ofile
    }
  },
  error = function(e) NULL
)
script_dir <- if (!is.null(this_script) && nzchar(this_script)) {
  dirname(normalizePath(this_script))
} else {
  getwd()
}

# ---- Locate VCFs ----------------------------------------------------------

vcf_files <- sort(list.files(
  script_dir,
  pattern = "\\.vcf\\.gz$",
  full.names = TRUE
))
stopifnot(length(vcf_files) > 0)

# Strip ".vcf.gz" and the "_INDELintersect" suffix for clean sample names.
sample_names <- sub(
  "_INDELintersect$",
  "",
  sub("\\.vcf\\.gz$", "", basename(vcf_files))
)

# ---- Build per-sample ID83 / ID89 / ID476 catalogs ------------------------
#
# mSigSpectra::vcfs_to_catalogs() does not yet support ID types, so run the
# per-sample pipeline (read_vcf -> split_vcf -> annotate_vcf -> vcf_to_catalog)
# and cbind the resulting one-column catalogs.

# Reference is hs37d5 (see VCF headers). Load the BSgenome object both for
# annotation and for the out-of-bounds POS check below.
bsg <- BSgenome::getBSgenome("BSgenome.Hsapiens.1000genomes.hs37d5")
chr_lens <- GenomeInfoDb::seqlengths(bsg)
margin <- 1000L

id_types <- c("ID83", "ID89", "ID476")
per_sample_cats <- setNames(vector("list", length(id_types)), id_types)
for (t in id_types) {
  per_sample_cats[[t]] <- vector("list", length(vcf_files))
}

for (i in seq_along(vcf_files)) {
  s <- sample_names[i]
  raw <- mSigSpectra::read_vcf(vcf_files[i], filter = "PASS", name_of_vcf = s)
  sp <- suppressWarnings(mSigSpectra::split_vcf(raw, name_of_vcf = s))
  id_vcf <- sp$ID

  # Drop variants whose POS lies outside the reference chromosome bounds with
  # enough margin for annotate_vcf()'s sequence-context fetch. At least one
  # input VCF contains a chr17 variant past the end of GRCh37/hs37d5 chr17.
  if (nrow(id_vcf) > 0L) {
    chroms <- as.character(id_vcf$CHROM)
    lens <- chr_lens[chroms]
    indel_w <- pmax(nchar(id_vcf$REF), nchar(id_vcf$ALT))
    ok <- !is.na(lens) &
      id_vcf$POS > margin &
      (id_vcf$POS + indel_w + margin) < lens
    dropped <- sum(!ok)
    if (dropped > 0L) {
      message(sprintf("%s: dropping %d out-of-bounds variant(s)", s, dropped))
      id_vcf <- id_vcf[ok, , drop = FALSE]
    }
  }

  ann <- mSigSpectra::annotate_id_vcf(
    id_vcf,
    ref_genome = bsg,
    name_of_vcf = s
  )

  for (t in id_types) {
    per_sample_cats[[t]][[
      i
    ]] <- suppressWarnings(mSigSpectra::vcf_to_id_catalog(
      ann,
      type = t,
      ref_genome = bsg,
      region = "genome",
      clip_le_9 = TRUE,
      FILTER_PASS = TRUE,
      sample_name = s
    ))
  }
}

id83 <- mSigSpectra::cbind_catalogs(per_sample_cats$ID83)
id89 <- mSigSpectra::cbind_catalogs(per_sample_cats$ID89)
id476 <- mSigSpectra::cbind_catalogs(per_sample_cats$ID476)

# ---- Plot (fig1.R layout: one letter-portrait page per sample) ------------

source(file.path(dirname(script_dir), "one_page_83_89_476.R"))

num_peaks <- 4
base_size <- 9.5
page_w <- 8.5 # letter portrait width  (in)
page_h <- 11 # letter portrait height (in)

out_pdf <- file.path(script_dir, "aa_vcfs_spectra2.pdf")

cairo_pdf(out_pdf, width = page_w, height = page_h, onefile = TRUE)

first_page <- TRUE
for (s in sample_names) {
  if (!first_page) {
    grid::grid.newpage()
  }
  first_page <- FALSE
  one_page_83_89_476(
    sig_83 = id83[, s, drop = FALSE],
    sig_89 = id89[, s, drop = FALSE],
    sig_476 = id476[, s, drop = FALSE],
    title_83 = paste0(s, " — ID83"),
    title_89 = paste0(s, " — ID89"),
    title_476 = paste0(s, " — ID476"),
    num_peak_labels = num_peaks,
    base_size = base_size,
    page_h = page_h
  )
}

grDevices::dev.off()

message("Wrote: ", out_pdf)
