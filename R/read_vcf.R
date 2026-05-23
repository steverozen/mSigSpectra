# Default FILTER values that are treated as "passing" by mSigSpectra.
# Union of the three ICAMS defaults: "PASS" (Strelka / Mutect), "." (Freebayes),
# and the empty string (missing value).
.default_pass_filter <- c("PASS", ".", "")

#' Read a VCF file into a data.table, caller-agnostically
#'
#' Reads the *body* of a VCF file (lines after `#CHROM`) into a `data.table`.
#' The resulting table has whatever columns the VCF has (`CHROM`, `POS`, `ID`,
#' `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`, and optionally `FORMAT` plus one
#' or more sample columns), with `#CHROM` renamed to `CHROM`.
#'
#' **Caller-agnostic.** `read_vcf()` does not know or care which variant
#' caller produced the VCF. It does not parse `FORMAT`/sample columns and
#' does not extract VAF or read depth. The only caller-dependent semantics
#' is the default value of the `filter` argument (see below).
#'
#' Uses [data.table::fread()] with
#' `check.names=FALSE`, `na.strings = ""``,
#' `fill = TRUE`, to parse the VCF body. Handles uncompressed
#' and gzipped files; does not handle bgzipped/tabix.
#'
#' @param file Path or URL to the VCF file.
#' @param filter Controls which rows are kept based on the `FILTER` column:
#'   * `TRUE` (default): keep rows where `FILTER %in% c("PASS", ".", "")` —
#'     the union of passing values across common callers.
#'   * `FALSE` or `NULL`: keep all rows.
#'   * A character vector: keep rows where `FILTER %in% filter`.
#'   Rows with no `FILTER` column are always kept; a warning is emitted
#'   when `filter` is non-trivial in that case.
#' @param name_of_vcf Optional name for the VCF, used only for warning /
#'   error messages. Defaults to the filename with extension stripped.
#'
#' @return A `data.table` with one row per variant. The name of the first
#'   column is 'CHROM', not '#CHROM'. Other column names
#'
#' @export
read_vcf <- function(file, filter = TRUE, name_of_vcf = NULL) {
  if (is.null(name_of_vcf)) {
    name_of_vcf <- tools::file_path_sans_ext(basename(file))
  }

  df <- read_vcf_fread(file)

  if (nrow(df) == 0L) {
    return(df)
  }

  filter_vals <- resolve_filter_arg(filter)
  if (!is.null(filter_vals)) {
    if (!"FILTER" %in% colnames(df)) {
      warning(
        "VCF ",
        dQuote(name_of_vcf),
        " has no FILTER column; returning all rows despite filter=",
        paste(filter_vals, collapse = ",")
      )
    } else {
      df <- df[df$FILTER %in% filter_vals, ]
    }
  }

  df
}

resolve_filter_arg <- function(filter) {
  if (isTRUE(filter)) {
    return(.default_pass_filter)
  }
  if (isFALSE(filter) || is.null(filter)) {
    return(NULL)
  }
  stopifnot(is.character(filter))
  filter
}

# Simple data.table::fread backend. Skips ## header lines by using
# skip = "#CHROM" to position fread at the column-header line.
read_vcf_fread <- function(file) {
  df1 <- tryCatch(
    suppressWarnings(data.table::fread(
      file,
      na.strings = "",
      skip = "#CHROM",
      fill = TRUE,
      check.names = FALSE,
    )),
    error = function(e) {
      stop(
        file,
        " does not appear to be a VCF file.\nDetails: ",
        conditionMessage(e)
      )
    }
  )

  if (nrow(df1) == 0L) {
    return(df1)
  }

  required <- c("#CHROM", "POS", "REF", "ALT")
  missing_cols <- setdiff(required, colnames(df1))
  if (length(missing_cols) > 0L) {
    stop(
      "Required VCF columns missing: ",
      paste(missing_cols, collapse = " ")
    )
  }

  data.table::setnames(df1, old = "#CHROM", new = "CHROM")
  df1$CHROM <- as.character(df1$CHROM)
  df1
}

# ---------------------------------------------------------------------------
# NOTE: vcfppR backend (removed for future development)
#
# An earlier version of read_vcf() supported `backend = "vcfppR"` to read
# bgzipped/tabix and remote VCF files via vcfppR::vcftable(). It was removed
# because vcfppR::vcftable() emitted warnings on some test VCFs and
# segfaulted on others (e.g. Strelka SBS), even with non-GT FORMAT tags.
#
# If a future maintainer wants to revisit this, the original backend body
# is preserved below for reference. To re-enable, restore `backend` as an
# argument to read_vcf(), re-add vcfppR to Suggests in DESCRIPTION, and
# diagnose the segfault in vcftable() (likely related to FORMAT parsing
# on files whose sample columns don't match what vcfppR expects).
#
# read_vcf_vcfppr <- function(file) {
#   if ("" == system.file(package = "vcfppR")) {
#     stop(
#       "\nPlease install vcfppR to use backend = 'vcfppR':\n",
#       "install.packages(\"vcfppR\")"
#     )
#   }
#   tbl <- vcfppR::vcftable(file)
#   dt <- data.table::data.table(
#     CHROM  = as.character(tbl$chr),
#     POS    = as.integer(tbl$pos),
#     ID     = tbl$id,
#     REF    = tbl$ref,
#     ALT    = tbl$alt,
#     QUAL   = tbl$qual,
#     FILTER = tbl$filter,
#     INFO   = tbl$info
#   )
#   dt
# }
# ---------------------------------------------------------------------------

#' Read multiple VCF files
#'
#' Thin batch wrapper around [read_vcf()].
#'
#' @param files Character vector of file paths / URLs.
#' @param ... Passed through to [read_vcf()].
#' @param names_of_vcfs Optional character vector of names (same length as
#'   `files`). Defaults to stripping extensions from basenames.
#'
#' @return A named list of `data.table`s, one per input file.
#'
#' @export
read_vcfs <- function(files, ..., names_of_vcfs = NULL) {
  if (is.null(names_of_vcfs)) {
    names_of_vcfs <- tools::file_path_sans_ext(basename(files))
  } else {
    stopifnot(length(names_of_vcfs) == length(files))
  }

  out <- Map(
    function(file, name) read_vcf(file, ..., name_of_vcf = name),
    files,
    names_of_vcfs
  )
  names(out) <- names_of_vcfs
  out
}
