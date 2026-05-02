#!/usr/bin/env Rscript
# Compare annotate_id_vcf() output against pre-existing annotated indel VCFs.
#
# For each of the N_FILES largest files in DATA_DIR:
#   1. Read the full old annotated VCF (plain CHROM header, no ##).
#   2. Strip to CHROM..FILTER (first 7 standard VCF cols).
#   3. Re-annotate with annotate_id_vcf().
#   4. Report columns present in old but absent in new, and vice-versa.
#   5. For each shared annotation column, report rows with differing values.

# ---- arguments --------------------------------------------------------------

p <- argparser::arg_parser("Compare annotate_id_vcf() output against pre-annotated indel VCFs")
p <- argparser::add_argument(p, "--n-files",  type = "integer", default = 1L,
                             help = "Number of largest files to process")
p <- argparser::add_argument(p, "--max-rows", type = "integer", default = NA_integer_,
                             help = "Cap on rows per file (omit for all rows)")
args <- argparser::parse_args(p)

DATA_DIR   <- "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs"
N_FILES    <- args$n_files
REF_GENOME <- "GRCh37"
MAX_ROWS   <- if (is.na(args$max_rows)) NULL else args$max_rows

# Columns expected to be absent from the new annotation — not reported as missing.
# Includes: VCF infrastructure not re-added, old dead-code columns, and
# per-sample columns (sample name varies per file; add them here as encountered).
OK_MISSING <- c(
  "INFO",
  "FORMAT",
  "CPCT02450014T",   # sample-name column (varies per file)
  "VAF",
  "read.depth",
  "seq.to.check",    # produced by old annotation pipeline, not by current code
  "prev_COSMIC_83"   # regression-test copy of old COSMIC_83 string
)

# ---- helpers ----------------------------------------------------------------

# Columns that are inputs or infrastructure — never compared.
.skip_always <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
)

compare_annotation <- function(old_vcf, new_vcf, discarded) {
  old_cols     <- colnames(old_vcf)
  new_cols     <- colnames(new_vcf)
  missing_cols <- setdiff(old_cols, new_cols)
  extra_cols   <- setdiff(new_cols, old_cols)

  unexpected_missing <- setdiff(missing_cols, OK_MISSING)
  cat("\nColumns in OLD but not in NEW — unexpected (",
      length(unexpected_missing), "):\n", sep = "")
  print(unexpected_missing)
  cat("\nColumns in OLD but not in NEW — expected/OK (",
      length(intersect(missing_cols, OK_MISSING)), "):\n", sep = "")
  print(intersect(missing_cols, OK_MISSING))
  cat("\nColumns in NEW but not in OLD (", length(extra_cols), "):\n", sep = "")
  print(extra_cols)

  n_discarded <- if (is.null(discarded)) 0L else nrow(discarded)
  cat(sprintf("\nRows discarded by justify_id_vcf: %d\n", n_discarded))

  # Match rows by CHROM + POS + REF + ALT.
  # The old VCFs are already justified (pos_shift = 0 typical), so POS matches.
  old_key <- do.call(paste, old_vcf[, c("CHROM", "POS", "REF", "ALT"), with = FALSE])
  new_key <- do.call(paste, new_vcf[, c("CHROM", "POS", "REF", "ALT"), with = FALSE])
  old_idx <- which(old_key %in% new_key)
  new_idx <- match(old_key[old_idx], new_key)

  cat(sprintf("Matched rows for value comparison: %d / %d\n\n",
              length(old_idx), nrow(old_vcf)))

  # Columns to compare: shared, excluding input/infrastructure, missing, extra,
  # and any sample-name column (identified as not being in the standard set).
  known_cols   <- union(.skip_always, union(missing_cols, extra_cols))
  compare_cols <- setdiff(intersect(old_cols, new_cols), known_cols)

  for (col in compare_cols) {
    old_vals <- old_vcf[[col]][old_idx]
    new_vals <- new_vcf[[col]][new_idx]
    diffs    <- which(
      !((old_vals == new_vals) | (is.na(old_vals) & is.na(new_vals)))
    )
    if (length(diffs) == 0L) {
      cat(sprintf("  OK   '%s': all %d rows match\n", col, length(old_idx)))
    } else {
      cat(sprintf("  DIFF '%s': %d / %d rows differ\n",
                  col, length(diffs), length(old_idx)))
      for (i in diffs[seq_len(min(3L, length(diffs)))]) {
        ri <- old_idx[i]
        cat(sprintf(
          "    CHROM=%s POS=%s REF=%s ALT=%s  old=<%s>  new=<%s>\n",
          old_vcf$CHROM[ri], old_vcf$POS[ri],
          old_vcf$REF[ri],   old_vcf$ALT[ri],
          as.character(old_vals[i]),
          as.character(new_vals[i])
        ))
      }
    }
  }
}

# ---- main -------------------------------------------------------------------

if (TRUE) {
  devtools::load_all()

  all_files <- list.files(
    path.expand(DATA_DIR),
    pattern    = "\\.annotated\\.indel\\.vcf\\.gz$",
    full.names = TRUE
  )
  sizes     <- file.info(all_files)$size
  top_files <- all_files[order(sizes, decreasing = TRUE)[seq_len(N_FILES)]]

  for (f in top_files) {
    cat(strrep("=", 72), "\n")
    cat("File:", basename(f), "\n")

    old_vcf <- data.table::fread(f, sep = "\t", quote = "")
    if (!is.null(MAX_ROWS)) {
      old_vcf <- old_vcf[seq_len(min(MAX_ROWS, nrow(old_vcf)))]
    }
    cat(sprintf("Rows read: %d   Cols: %d\n", nrow(old_vcf), ncol(old_vcf)))

    filter_idx <- which(colnames(old_vcf) == "FILTER")
    stopifnot(length(filter_idx) == 1L)
    stripped <- old_vcf[, seq_len(filter_idx), with = FALSE]

    result    <- suppressWarnings(
      annotate_id_vcf(stripped, ref_genome = REF_GENOME, explain_indels = 0)
    )

    compare_annotation(old_vcf, result$annotated.vcf, result$discarded.variants)
  }
}
