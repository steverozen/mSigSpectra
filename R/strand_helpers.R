#' Reverse complement a character vector of DNA strings
#'
#' Thin vectorized wrapper around
#' [Biostrings::reverseComplement()][Biostrings::reverse]. Handles IUPAC
#' ambiguity codes.
#'
#' @param x A character vector of DNA sequences.
#'
#' @return A character vector of reverse complements.
#'
#' @export
revc <- function(x) {
  if (length(x) == 0L) return(character(0L))
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

#' Normalize SBS1536 pentanucleotide + ALT strings to pyrimidine form
#'
#' Input strings are 6 characters: a 5-base pentanucleotide context followed
#' by a 1-base ALT (e.g. `"ATGCTT"` = `ATGCT>T`). If the center of the
#' pentanucleotide (position 3) is `A` or `G`, the pentanucleotide and the
#' ALT are reverse-complemented so that the center becomes `C` or `T` —
#' the pyrimidine form used in canonical SBS96/1536 row names.
#'
#' @param x A character vector of 6-letter strings.
#'
#' @keywords internal
pyr_penta <- function(x) {
  stopifnot(all(nchar(x) == 6L))
  center <- substr(x, 3L, 3L)
  needs_rc <- center %in% c("A", "G")
  out <- x
  if (any(needs_rc)) {
    out[needs_rc] <- paste0(
      revc(substr(x[needs_rc], 1L, 5L)),
      revc(substr(x[needs_rc], 6L, 6L))
    )
  }
  out
}

#' Reverse-complement stranded SBS96 4-character class strings
#'
#' Input is a 4-character string where characters 1-3 are the trinucleotide
#' context and character 4 is the ALT (e.g. `"AATC"` = `AAT>ACT`).
#'
#' @param x A character vector of 4-letter strings.
#'
#' @keywords internal
revc_sbs96 <- function(x) {
  stopifnot(all(nchar(x) == 4L))
  paste0(revc(substr(x, 1L, 3L)), revc(substr(x, 4L, 4L)))
}

#' Reverse-complement stranded DBS144 4-character class strings
#'
#' Input is a 4-character string where characters 1-2 are the REF dinucleotide
#' and characters 3-4 are the ALT (e.g. `"AATC"` = `AA>TC`).
#'
#' @param x A character vector of 4-letter strings.
#'
#' @keywords internal
revc_dbs144 <- function(x) {
  stopifnot(all(nchar(x) == 4L))
  paste0(revc(substr(x, 1L, 2L)), revc(substr(x, 3L, 4L)))
}
