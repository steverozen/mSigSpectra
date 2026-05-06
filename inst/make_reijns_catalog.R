# make_reijns_catalog.R
#
# Read per-sample Reijns RPE1 VCF files, annotate their indels with
# mSigSpectra (hg38), and write two Koh-476 indel catalogs:
#   - keep_high_R: all indels retained (clip_le_9 = FALSE)
#   - drop_high_R: R > 9 single-base T/C indels clipped (clip_le_9 = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(mSigSpectra)
  library(BSgenome)
  library(GenomeInfoDb)
})

#' Build Koh-476 indel catalogs from Reijns RPE1 per-sample VCFs
#'
#' Reads every \code{reijns_cell_*.vcf} file in \code{vcf_dir}, annotates
#' indels with \code{\link{annotate_id_vcf}}, and calls
#' \code{\link{vcf_to_id_catalog}} twice — once with \code{clip_le_9 = FALSE}
#' (retain R > 9 variants) and once with \code{clip_le_9 = TRUE} (drop them) —
#' producing two tab-separated catalog files with one column per sample.
#'
#' @param vcf_dir Directory containing \code{reijns_cell_*.vcf} files.
#' @param ref_genome BSgenome object or alias; \code{"hg38"} by default.
#' @param out_keep Output TSV path for the catalog with R > 9 kept.
#' @param out_drop Output TSV path for the catalog with R > 9 dropped.
make_reijns_catalog <- function(
  vcf_dir = NULL,
  ref_genome = "hg38",
  out_keep = NULL,
  out_drop = NULL
) {
  if (is.null(vcf_dir)) {
    ca <- commandArgs(trailingOnly = FALSE)
    m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
    vcf_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()
  }
  if (is.null(out_keep)) {
    out_keep <- file.path(vcf_dir, "reijns_cells_catalog_keep_high_R.tsv")
  }
  if (is.null(out_drop)) {
    out_drop <- file.path(vcf_dir, "reijns_cells_catalog_drop_high_R.tsv")
  }

  vcf_files <- sort(list.files(
    vcf_dir,
    pattern = "reijns_cell.*\\.vcf$",
    full.names = TRUE
  ))
  if (length(vcf_files) == 0L) {
    stop("No reijns_cell*.vcf files found in: ", vcf_dir)
  }
  message("Found ", length(vcf_files), " VCF file(s) in ", vcf_dir)

  bsg <- mSigSpectra:::normalize_genome_arg(ref_genome)
  chr_lens <- GenomeInfoDb::seqlengths(bsg)
  margin <- 1000L

  cols_keep <- vector("list", length(vcf_files))
  cols_drop <- vector("list", length(vcf_files))

  for (i in seq_along(vcf_files)) {
    f <- vcf_files[i]
    s <- sub("reijns_cell_", "", sub("\\.vcf$", "", basename(f)))
    message("Processing ", basename(f))

    raw <- mSigSpectra::read_vcf(f, filter = "PASS", name_of_vcf = s)
    sp <- suppressWarnings(mSigSpectra::split_vcf(raw, name_of_vcf = s))
    id_vcf <- sp$ID

    if (nrow(id_vcf) == 0L) {
      message("  no indels — skipping")
      next
    }

    chroms <- as.character(id_vcf$CHROM)
    lens <- chr_lens[chroms]
    indel_w <- pmax(nchar(id_vcf$REF), nchar(id_vcf$ALT))
    ok <- !is.na(lens) &
      id_vcf$POS > margin &
      (id_vcf$POS + indel_w + margin) < lens
    dropped <- sum(!ok)
    if (dropped > 0L) {
      message("  dropping ", dropped, " out-of-bounds variant(s)")
      id_vcf <- id_vcf[ok, , drop = FALSE]
    }
    if (nrow(id_vcf) == 0L) {
      next
    }

    ann <- mSigSpectra::annotate_id_vcf(
      id_vcf,
      ref_genome = bsg,
      name_of_vcf = s
    )

    cols_keep[[i]] <- mSigSpectra::vcf_to_id_catalog(
      ann,
      type        = "ID476",
      ref_genome  = bsg,
      sample_name = s,
      FILTER_PASS = FALSE,
      clip_le_9   = FALSE
    )

    cols_drop[[i]] <- mSigSpectra::vcf_to_id_catalog(
      ann,
      type        = "ID476",
      ref_genome  = bsg,
      sample_name = s,
      FILTER_PASS = FALSE,
      clip_le_9   = TRUE
    )
  }

  write_id476_tsv(cols_keep, out_keep)
  write_id476_tsv(cols_drop, out_drop)

  invisible(list(keep = out_keep, drop = out_drop))
}

write_id476_tsv <- function(cols, path) {
  valid <- !vapply(cols, is.null, logical(1))
  if (!any(valid)) {
    warning("No samples with indels; catalog not written to ", path)
    return(invisible(NULL))
  }
  mat <- do.call(cbind, cols[valid])
  dt <- data.table::data.table(
    MutationType = rownames(mat),
    data.table::as.data.table(as.data.frame(mat))
  )
  data.table::fwrite(dt, path, sep = "\t")
  message(
    "Wrote ", path,
    " (", nrow(mat), " mutation types × ", sum(valid), " samples)"
  )
}
