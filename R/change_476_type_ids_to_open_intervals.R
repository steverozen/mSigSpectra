#' Change 476-type indel category identifiers to use right-open repeat intervals
#'
#' Replaces bounded repeat-count suffixes like `:R(X,9)` with right-open
#' equivalents like `:R(X,)`, making the classification agnostic as to whether
#' repeats longer than 9 were discarded from the input data.
#'
#' @param type_476_indel_type_identifiers Character vector of 476-type indel
#'   category identifiers, e.g. as returned by
#'   [categorize_indels_in_vcf()].
#'
#' @return Character vector the same length as the input, with `:R(X,9)` at the
#'   end of each string replaced by `:R(X,)`.
#'
#' @examples
#' change_476_type_ids_to_open_intervals(
#'   c("Del(C):Ins(C):R(5,9)", "Del(T):R(3,5)", "Ins(C):R(5,9)")
#' )
#'
#' @export
change_476_type_ids_to_open_intervals = function(
  type_476_indel_type_identifiers
) {
  gsub(":R\\((.),9\\)$", ":R(\\1,)", x = type_476_indel_type_identifiers)
}
