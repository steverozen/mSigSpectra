#' Normalize SBS1536 pentanucleotide + ALT strings to pyrimidine form
#'
#' Input strings are 6 characters: a 5-base pentanucleotide context followed
#' by a 1-base ALT (e.g. `"ATGCTT"` = `ATGCT>T`). If the center of the
#' pentanucleotide (position 3) is `A` or `G`, the pentanucleotide and the
#' ALT are reverse-complemented so that the center becomes `C` or `T` —
#' the pyrimidine form used in canonical SBS96/1536 row names.
#'
#' @param mutstring A character vector of 6-letter strings.
#'
#' @keywords internal
pyr_penta <- function(mutstring) {
  stopifnot(all(nchar(mutstring) == 6L))
  center <- substr(mutstring, 3L, 3L)
  needs_rc <- center %in% c("A", "G")
  out <- mutstring
  if (any(needs_rc)) {
    out[needs_rc] <- paste0(
      fastrc::fast_rc(substr(mutstring[needs_rc], 1L, 5L)),
      fastrc::fast_rc(substr(mutstring[needs_rc], 6L, 6L))
    )
  }
  out
}

#' Reverse complement strings that represent stranded SBSs
#'
#' Input is a 4-character string where characters 1-3 are the trinucleotide
#' context and character 4 is the ALT (e.g. `"AATC"` = `AAT>ACT`).
#' Returns the reverse complement of the first 3 characters concatenated
#' with the reverse complement of the last character,
#' e.g. `"AATC"` returns `"ATTG"`.
#'
#' @param mutstring A character vector of 4-letter strings.
#'
#' @keywords internal
revc_sbs96 <- function(mutstring) {
  stopifnot(all(nchar(mutstring) == 4L))
  paste0(fastrc::fast_rc(substr(mutstring, 1L, 3L)), fastrc::fast_rc(substr(mutstring, 4L, 4L)))
}

#' Reverse complement strings that represent stranded DBSs
#'
#' Input is a 4-character string where characters 1-2 are the REF dinucleotide
#' and characters 3-4 are the ALT (e.g. `"AATC"` = `AA>TC`).
#' Returns the reverse complement of the first 2 characters concatenated
#' with the reverse complement of the last 2 characters,
#' e.g. `"AATC"` returns `"TTGA"`.
#'
#' @param mutstring A character vector of 4-letter strings.
#'
#' @keywords internal
revc_dbs144 <- function(mutstring) {
  stopifnot(all(nchar(mutstring) == 4L))
  paste0(fastrc::fast_rc(substr(mutstring, 1L, 2L)), fastrc::fast_rc(substr(mutstring, 3L, 4L)))
}
