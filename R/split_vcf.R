#' Split a mixed-mutation VCF into SBS / DBS / ID sub-tables
#'
#' Classifies each row by `REF`/`ALT` length alone:
#' * SBS (single base substitution): `nchar(REF) == 1 && nchar(ALT) == 1`
#' * DBS (double base substitution): `nchar(REF) == 2 && nchar(ALT) == 2`
#' * ID (indel):                     `nchar(REF) != nchar(ALT)`
#' * Other (e.g. 3+bp substitutions): discarded with a recorded reason.
#'
#' This is **caller-agnostic**. In particular, mSigSpectra does *not*
#' merge adjacent SBSs into DBSs based on VAF similarity (which ICAMS did
#' via `SplitOneVCF` for Strelka-style VCFs). Users who want that behavior
#' should apply it as a post-processing step with their own VAF column.
#'
#' @param vcf A VCF as a data.frame / data.table with at least `REF` and
#'   `ALT` columns.
#' @param name_of_vcf Optional VCF name used in warning / error messages.
#'
#' @return A list with elements
#'   * `SBS`: `data.table` of SBS rows.
#'   * `DBS`: `data.table` of DBS rows.
#'   * `ID`: `data.table` of indel rows.
#'   * `discarded`: `data.table` of rows that did not fit any of the above,
#'     with a `discarded.reason` column — **`NULL` when no rows were
#'     discarded**.
#' @md
#'
#' @export
split_vcf <- function(vcf, name_of_vcf = NULL) {
  vcf <- data.table::as.data.table(vcf)
  empty <- vcf[0, ]

  if (nrow(vcf) == 0L) {
    return(list(SBS = empty, DBS = empty, ID = empty, discarded = NULL))
  }

  ref_n <- nchar(vcf$REF)
  alt_n <- nchar(vcf$ALT)

  sbs_i <- ref_n == 1L & alt_n == 1L
  dbs_i <- ref_n == 2L & alt_n == 2L
  id_i  <- ref_n != alt_n
  other_i <- !(sbs_i | dbs_i | id_i)

  discarded <- if (any(other_i)) {
    d <- vcf[other_i, ]
    d$discarded.reason <-
      "REF/ALT did not match SBS (1/1), DBS (2/2), or ID (unequal length)"
    warning(
      "In VCF ",
      ifelse(is.null(name_of_vcf), "", dQuote(name_of_vcf)),
      " ", nrow(d), " row out of ", nrow(vcf),
      " had REF/ALT combinations that are not SBS/DBS/ID and were discarded. ",
      "See split_vcf()$discarded for details."
    )
    d
  } else {
    NULL
  }

  list(
    SBS = vcf[sbs_i, ],
    DBS = vcf[dbs_i, ],
    ID  = vcf[id_i, ],
    discarded = discarded
  )
}
