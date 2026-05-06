# make_reijns_rosetta_csv.R
#
# Read per-sample Reijns RPE1 VCF files from a directory, annotate their
# indels with mSigSpectra (hg38), and build a flat rosetta-stone CSV.
#
# Adapted from make_rosetta_csv.R; the only structural difference is that
# it reads plain .vcf files (not .vcf.gz) produced by make_reijns_sample_vcfs.R.
#
# With open_interval_format = TRUE in gen_476_type_string.R, single-base T/C
# repeats with R >= 9 all produce "R(9,)" in Koh_476, so R > 9 rows land in
# the R(9,) catalog bin rather than being "uncapped".  The drop_high_R option
# therefore acts on the numeric R column, not on the Koh_476 string.

suppressPackageStartupMessages({
  library(data.table)
  library(mSigSpectra)
  library(BSgenome)
  library(GenomeInfoDb)
})

# ---- helpers --------------------------------

shorten_long_visual <- function(x, flank_5 = 5L, flank_3 = 20L) {
  out <- x
  ok <- !is.na(x)
  parts <- strsplit(x[ok], " ", fixed = TRUE)
  ok_three <- vapply(parts, length, integer(1)) == 3L
  idx <- which(ok)[ok_three]
  good <- parts[ok_three]
  pre <- vapply(good, `[`, character(1), 1)
  mid <- vapply(good, `[`, character(1), 2)
  post <- vapply(good, `[`, character(1), 3)
  out[idx] <- paste(
    substr(pre, pmax(1L, nchar(pre) - flank_5 + 1L), nchar(pre)),
    mid,
    substr(post, 1L, flank_3)
  )
  out
}

is_single_tc_class <- function(x) {
  grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R[^]]+\\][ACGT]$", x)
}

# ---- main function ------------------------------------------------------------

#' Build a rosetta-stone CSV from Reijns RPE1 per-sample indel VCFs
#'
#' Reads every \code{reijns_cell_*.vcf} file in \code{vcf_dir}, annotates
#' indels with \code{\link{annotate_id_vcf}} using the hg38 genome, then builds
#' a flat CSV mapping 476-type / 89-type / 83-type identifiers with example
#' \code{long_visual} strings per pair.  The per-sample VCF files are
#' produced by \code{make_reijns_sample_vcfs.R}.
#'
#' @param vcf_dir Directory containing \code{reijns_cell_*.vcf} files.
#' @param ref_genome BSgenome object or alias; \code{"hg38"} by default.
#' @param out_path Output CSV path.
#' @param drop_high_R If \code{TRUE}, drop rows where the numeric repeat count
#'   R \eqn{>} 9 for single-base T/C classes.  Either way the count is reported.
#' @param flank_5 Characters to keep from the 5' flank in \code{long_visual}.
#' @param flank_3 Characters to keep from the 3' flank in \code{long_visual}.
#' @param n_examples Maximum example rows per (476-type, 89-type) pair for
#'   non-single-base T/C classes.
#' @param n_singletc Example rows per pair for single-base T/C classes.
make_reijns_rosetta_csv <- function(
  vcf_dir = NULL,
  ref_genome = "hg38",
  out_path = NULL,
  drop_high_R = FALSE,
  flank_5 = 5L,
  flank_3 = 20L,
  n_examples = 20L,
  n_singletc = 1L
) {
  if (is.null(vcf_dir)) {
    ca <- commandArgs(trailingOnly = FALSE)
    m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
    vcf_dir <- if (length(m)) dirname(normalizePath(m[1])) else getwd()
  }
  if (is.null(out_path)) {
    out_path <- file.path(vcf_dir, "reijns_rosetta.csv")
  }

  keep_cols <- c(
    "Koh_476",
    "Koh_89",
    "COSMIC_83",
    "long_visual",
    "ins_or_del_seq",
    "U_seq",
    "U_seq_count_in_indel_seq",
    "R",
    "mh",
    "unit",
    "unit_length",
    "internal_rep",
    "internal_reps",
    "spacer",
    "spacer_length",
    "prime3_rep",
    "prime3_reps",
    "original_reps"
  )

  # ---- 1. Read and annotate VCFs -------------------------------------------

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

  all_rows <- vector("list", length(vcf_files))

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
    )$annotated.vcf

    missing <- setdiff(keep_cols, colnames(ann))
    for (col in missing) {
      ann[[col]] <- NA_character_
    }
    ann <- as.data.table(ann)[, keep_cols, with = FALSE]

    ann[, long_visual := shorten_long_visual(long_visual, flank_5, flank_3)]
    all_rows[[i]] <- ann
  }

  dt <- data.table::rbindlist(all_rows, fill = TRUE)
  message("Total rows after annotation: ", nrow(dt))

  # Drop rows whose Koh_476 is not a valid catalog entry (e.g. very long
  # homopolymers where the context window is too narrow cause the regex to
  # backtrack and produce an impossible right-flank base).
  valid_476 <- mSigSpectra::catalog_row_order()$ID476
  bad_476 <- !(dt$Koh_476 %in% valid_476)
  if (any(bad_476)) {
    message(
      sum(bad_476),
      " row(s) with unrecognised Koh_476 dropped: ",
      paste(unique(dt$Koh_476[bad_476]), collapse = ", ")
    )
    dt <- dt[!bad_476]
  }

  # ---- 2. R > 9 report / filter (uses numeric R column) --------------------
  # With open_interval_format = TRUE, single-base T/C R >= 9 all land in the
  # R(9,) catalog bin, so Koh_476 strings no longer distinguish R = 9 from
  # R = 13.  Use the numeric R column to count and optionally remove them.

  is_high_R <- !is.na(dt$R) &
    as.integer(dt$R) > 9L &
    is_single_tc_class(dt$Koh_476)
  n_high_R <- sum(is_high_R)
  if (drop_high_R) {
    message(n_high_R, " rows with R > 9 removed")
    dt <- dt[!is_high_R]
  } else {
    message(n_high_R, " rows with R > 9 retained")
  }

  # ---- 3. Drop strand-flipped single-base T/C rows (non-canonical base) -----

  is_single_tc <- is_single_tc_class(dt$Koh_476)
  drop_flip <- is_single_tc & !grepl(" <[CT]>", dt$long_visual)
  message(
    sum(drop_flip),
    " rows dropped: strand-flipped single-base T/C (non-canonical base hidden)"
  )
  dt <- dt[!drop_flip]

  # ---- 4. Deduplicate on classification + example string --------------------

  dt_unique <- unique(dt, by = c("Koh_476", "Koh_89", "long_visual"))
  setorder(dt_unique, Koh_476, Koh_89)

  # ---- 5. Derive all unique (Koh_476, Koh_89) pairs -------------------------

  pairs <- unique(dt_unique[, .(Koh_476, Koh_89)])

  # ---- 6. Select up to n_examples per pair (round-robin by ins_or_del_seq) --

  extra_cols <- c(
    "ins_or_del_seq",
    "U_seq",
    "U_seq_count_in_indel_seq",
    "R",
    "mh",
    "unit",
    "unit_length",
    "internal_rep",
    "internal_reps",
    "spacer",
    "spacer_length",
    "prime3_rep",
    "prime3_reps",
    "original_reps"
  )

  examples <- dt_unique[,
    {
      if (is_single_tc_class(.BY$Koh_476)) {
        head(.SD, n_singletc)
      } else {
        seqs <- unique(ins_or_del_seq)
        if (length(seqs) >= n_examples) {
          .SD[match(seqs[seq_len(n_examples)], ins_or_del_seq)]
        } else {
          seq_rank <- match(ins_or_del_seq, seqs)
          within_rank <- ave(seq_rank, seq_rank, FUN = seq_along)
          head(.SD[order(within_rank, seq_rank)], n_examples)
        }
      }
    },
    by = .(Koh_476, Koh_89),
    .SDcols = c("COSMIC_83", "long_visual", extra_cols)
  ]

  # ---- 7. Aggregate COSMIC_83 as "; "-separated string per pair -------------

  cosmic_per_pair <- dt_unique[,
    .(COSMIC_83 = paste(sort(unique(COSMIC_83)), collapse = "; ")),
    by = .(Koh_476, Koh_89)
  ]

  # ---- 8. Merge pairs + examples + aggregated COSMIC_83 ---------------------

  doc <- merge(
    pairs,
    examples[, c("Koh_476", "Koh_89", "long_visual", extra_cols), with = FALSE],
    by = c("Koh_476", "Koh_89"),
    all.x = TRUE,
    sort = FALSE
  )
  doc <- merge(
    doc,
    cosmic_per_pair,
    by = c("Koh_476", "Koh_89"),
    all.x = TRUE,
    sort = FALSE
  )
  doc[, example_n := seq_len(.N), by = .(Koh_476, Koh_89)]

  # ---- 9. Collapse spaces in long_visual ------------------------------------

  doc[, long_visual := gsub(" ", "", long_visual, fixed = TRUE)]

  # ---- 10. Order by canonical ID476 row order --------------------------------

  ord <- mSigSpectra::catalog_row_order()$ID476
  doc[, .row_ord := match(Koh_476, ord)]
  not_found <- unique(doc$Koh_476[is.na(doc$.row_ord)])
  if (length(not_found) > 0L) {
    warning(
      "Koh_476 values not in catalog_row_order()$ID476: ",
      paste(not_found, collapse = ", ")
    )
  }
  setorder(doc, .row_ord, Koh_89, example_n, na.last = TRUE)
  doc[, .row_ord := NULL]

  # ---- 11. Mark 476-type values that map to > 1 89-type with "*" ------------

  multi_476 <- doc[, data.table::uniqueN(Koh_89), by = Koh_476][
    V1 > 1L,
    Koh_476
  ]
  if (length(multi_476) > 0L) {
    message(
      "476-type classes mapping to multiple 89-types (marked with *): ",
      paste(multi_476, collapse = ", ")
    )
    doc[Koh_476 %in% multi_476, Koh_476 := paste0(Koh_476, "*")]
  }

  # ---- 12. Final column order and names -------------------------------------

  doc <- doc[,
    c(
      "Koh_476",
      "Koh_89",
      "COSMIC_83",
      "example_n",
      "long_visual",
      extra_cols
    ),
    with = FALSE
  ]

  setnames(
    doc,
    old = c("Koh_476", "Koh_89", "COSMIC_83", "example_n", "long_visual"),
    new = c("476-type", "89-type", "83-type", "example_n", "indel_in_context")
  )

  # ---- 13. Write CSV --------------------------------------------------------

  n_pairs <- data.table::uniqueN(doc, by = c("476-type", "89-type"))
  data.table::fwrite(doc, out_path)
  message(
    "Wrote ",
    out_path,
    " (",
    nrow(doc),
    " rows, ",
    n_pairs,
    " (476-type, 89-type) pairs)"
  )

  invisible(doc)
}
