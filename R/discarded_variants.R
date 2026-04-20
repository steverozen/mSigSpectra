#' Remove rows that look like stray "#CHROM" header repeats
#'
#' Occasionally a VCF body contains a spurious row whose `CHROM` column
#' equals `"#CHROM"` (e.g. when a concatenated VCF keeps the header line
#' of each input). Drop those rows and record them in `discarded.variants`.
#'
#' @keywords internal
remove_rows_with_pound_sign <- function(df, name_of_vcf = NULL) {
  idx <- which(df$CHROM == "#CHROM")
  if (length(idx) == 0L) return(list(df = df))

  warning(
    "In VCF ",
    ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)), " ",
    length(idx), " row out of ", nrow(df),
    " had value '#CHROM' in column CHROM and were removed. ",
    "See discarded.variants in the return value for more details"
  )
  to_remove <- df[idx, ]
  to_remove$discarded.reason <- 'Chromosome name is "#CHROM"'
  list(df = df[-idx, ], discarded.variants = to_remove)
}

#' Remove variants with duplicated CHROM+POS
#'
#' Two passes:
#' 1. Rows with identical (CHROM, POS, REF, ALT): keep one copy, discard the rest.
#' 2. Rows sharing (CHROM, POS, REF) but different ALT: discard all (treated as
#'    unresolved multiallelic / inconsistent records).
#'
#' @keywords internal
remove_rows_with_duplicated_chrom_and_pos <- function(df, name_of_vcf = NULL) {
  discarded <- df[0, ]

  dups_exact <- which(duplicated(df[, c("CHROM", "POS", "REF", "ALT")]))
  if (length(dups_exact) > 0L) {
    warning(
      "In VCF ",
      ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)), " ",
      2L * length(dups_exact), " row out of ", nrow(df),
      " had same CHROM, POS, REF, ALT and only one copy is kept. ",
      "See discarded.variants in the return value for more details"
    )
    to_remove <- df[dups_exact, ]
    to_remove$discarded.reason <-
      "Variant with same CHROM, POS, REF and ALT as another variant"
    discarded <- rbind(discarded, to_remove)
    df <- df[-dups_exact, ]
  }

  dups_alt <- which(duplicated(df[, c("CHROM", "POS", "REF")]))
  if (length(dups_alt) > 0L) {
    warning(
      "In VCF ",
      ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)), " ",
      2L * length(dups_alt), " row out of ", nrow(df),
      " had same CHROM, POS, REF but different ALT and were removed. ",
      "See discarded.variants in the return value for more details"
    )
    dups_fwd <- which(duplicated(df[, c("CHROM", "POS", "REF")], fromLast = TRUE))
    to_remove <- df[c(dups_alt, dups_fwd), ]
    to_remove$discarded.reason <-
      "Variant with same CHROM, POS, REF but different ALT"
    discarded <- rbind(discarded, to_remove)
    df <- df[-c(dups_alt, dups_fwd), ]
  }

  if (nrow(discarded) == 0L) list(df = df)
  else list(df = df, discarded.variants = discarded)
}

#' Check a VCF for common variant-level problems and remove the offenders
#'
#' Removes:
#' * Rows with identical REF and ALT.
#' * Stray `#CHROM` header repeats.
#' * Duplicated (CHROM, POS, REF, ALT) rows (keeping one copy).
#' * Multiple-ALT rows (comma-separated ALT field).
#' * Non-standard chromosome names (or those outside
#'   `chr_names_to_process` when supplied).
#' * Substitutions of length > 2 (e.g. `ACT>TGA`).
#' * Complex indels (REF[1] != ALT[1]).
#' * Wrong DBS rows where REF and ALT share a base at the same position.
#' * Variants whose REF base is not in `{A, C, G, T}`.
#'
#' Each discarded row gains a `discarded.reason` column. Returns a list
#' with `df` (retained rows) and optionally `discarded.variants`
#' (discarded rows).
#'
#' @param vcf A VCF as a data.frame / data.table.
#' @param name_of_vcf Optional name, used in warning messages.
#' @param chr_names_to_process Optional character vector of chromosome
#'   names to keep (overrides the default non-standard-contig filter).
#'
#' @export
check_and_remove_discarded_variants <- function(vcf,
                                                name_of_vcf = NULL,
                                                chr_names_to_process = NULL) {
  if (nrow(vcf) == 0L) return(list(df = vcf))

  discarded <- vcf[0, ]

  # Same REF and ALT
  idx <- which(vcf$REF == vcf$ALT)
  if (length(idx) > 0L) {
    to_remove <- vcf[idx, ]
    to_remove$discarded.reason <- "Variant with same REF and ALT"
    discarded <- rbind(discarded, to_remove)
    vcf <- vcf[-idx, ]
  }

  ret <- remove_rows_with_pound_sign(vcf, name_of_vcf = name_of_vcf)
  vcf <- ret$df
  if (!is.null(ret$discarded.variants)) {
    discarded <- rbind(discarded, ret$discarded.variants)
  }

  ret <- remove_rows_with_duplicated_chrom_and_pos(vcf, name_of_vcf = name_of_vcf)
  vcf <- ret$df
  if (!is.null(ret$discarded.variants)) {
    discarded <- rbind(discarded, ret$discarded.variants)
  }

  if (is.null(chr_names_to_process)) {
    ret <- standard_chrom_name_new(vcf, name.of.VCF = name_of_vcf)
  } else {
    ret <- select_variants_by_chrom_name(
      vcf, chr.names.to.process = chr_names_to_process,
      name.of.VCF = name_of_vcf
    )
  }
  vcf <- ret$df
  if (!is.null(ret$discarded.variants)) {
    discarded <- rbind(discarded, ret$discarded.variants)
  }

  # Multiple-ALT rows
  multi_alt <- grep(",", vcf$ALT, fixed = TRUE)
  if (length(multi_alt) > 0L) {
    warning(
      "VCF ", ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)),
      " has variants with multiple alternative alleles which were discarded. ",
      "See discarded.variants in the return value for more details."
    )
    to_remove <- vcf[multi_alt, ]
    to_remove$discarded.reason <- "Variant with multiple alternative alleles"
    discarded <- rbind(discarded, to_remove)
    vcf <- vcf[-multi_alt, ]
  }

  # Substitutions of length > 2 (both REF and ALT > 2 and equal)
  long_sub <- which(nchar(vcf$REF) > 2L & nchar(vcf$ALT) == nchar(vcf$REF))
  if (length(long_sub) > 0L) {
    warning(
      "VCF ", ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)),
      " has variants involving three or more nucleotides which were discarded."
    )
    to_remove <- vcf[long_sub, ]
    to_remove$discarded.reason <- "Variant involves three or more nucleotides"
    discarded <- rbind(discarded, to_remove)
    vcf <- vcf[-long_sub, ]
  }

  # Complex indels (REF[1] != ALT[1] with unequal lengths)
  complex_idx <- which(
    nchar(vcf$REF) > 0L & nchar(vcf$ALT) > 0L &
      nchar(vcf$REF) != nchar(vcf$ALT) &
      substr(vcf$REF, 1L, 1L) != substr(vcf$ALT, 1L, 1L)
  )
  if (length(complex_idx) > 0L) {
    warning(
      "VCF ", ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)),
      " has complex indels which were discarded."
    )
    to_remove <- vcf[complex_idx, ]
    to_remove$discarded.reason <- "Complex indel"
    discarded <- rbind(discarded, to_remove)
    vcf <- vcf[-complex_idx, ]
  }

  # Wrong DBS rows: REF and ALT share a base at the same position
  wrong_dbs <- which(
    nchar(vcf$REF) == 2L & nchar(vcf$ALT) == 2L &
      (substr(vcf$REF, 1L, 1L) == substr(vcf$ALT, 1L, 1L) |
       substr(vcf$REF, 2L, 2L) == substr(vcf$ALT, 2L, 2L))
  )
  if (length(wrong_dbs) > 0L) {
    warning(
      "VCF ", ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)),
      " has wrong DBS variants which were discarded."
    )
    to_remove <- vcf[wrong_dbs, ]
    to_remove$discarded.reason <- "Wrong DBS variant"
    discarded <- rbind(discarded, to_remove)
    vcf <- vcf[-wrong_dbs, ]
  }

  # Ambiguous REF (not A/C/G/T at position 1)
  amb_ref <- which(!substr(vcf$REF, 1L, 1L) %in% c("A", "C", "G", "T"))
  if (length(amb_ref) > 0L) {
    warning(
      "VCF ", ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)),
      " has variants with ambiguous REF bases which were discarded."
    )
    to_remove <- vcf[amb_ref, ]
    to_remove$discarded.reason <- "Variant has ambiguous REF base"
    discarded <- rbind(discarded, to_remove)
    vcf <- vcf[-amb_ref, ]
  }

  if (nrow(discarded) == 0L) list(df = vcf)
  else list(df = vcf, discarded.variants = discarded)
}
