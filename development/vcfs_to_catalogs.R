#' Convenience wrapper: read, annotate, and catalog a set of VCF files
#'
#' Thin composition of [read_vcf()] + [split_vcf()] +
#' [annotate_sbs_or_dbs_vcf()] + [vcf_to_sbs_catalog()] for a list of
#' files.
#'
#' NOTE: currently only SBS catalog types are wired through; DBS and ID
#' types will raise an error. This file lives under `development/` and is
#' excluded from the built package via `.Rbuildignore`.
#'
#' @param files Character vector of VCF paths / URLs.
#' @param types Character vector of catalog types to produce (e.g.
#'   `c("SBS96", "SBS192")`). Each file contributes one column per type.
#' @param ref_genome A BSgenome object or character alias.
#' @param region One of `"genome"`, `"exome"`, `"transcript"`,
#'   `"unknown"`.
#' @param trans_ranges Optional transcript ranges (see
#'   [add_transcript_strand()]).
#' @param names_of_vcfs Optional character vector of names (same length as
#'   `files`); defaults to basenames without extension.
#'
#' @return A named list keyed by `types`. Each element is a multi-sample
#'   catalog matrix (columns = input files).
#'
#' @export
vcfs_to_catalogs <- function(files, types, ref_genome,
                             region = "genome",
                             trans_ranges = NULL,
                             names_of_vcfs = NULL) {
  if (is.null(names_of_vcfs)) {
    names_of_vcfs <- tools::file_path_sans_ext(basename(files))
  }
  stopifnot(length(names_of_vcfs) == length(files))

  needs_sbs <- any(types %in% c("SBS96", "SBS192", "SBS1536"))
  needs_dbs <- any(types %in% c("DBS78", "DBS136", "DBS144"))
  needs_id  <- any(types %in% c("ID83", "ID89", "ID166", "ID476"))

  per_sample <- lapply(seq_along(files), function(i) {
    raw <- read_vcf(files[i], name_of_vcf = names_of_vcfs[i])
    sp <- suppressWarnings(split_vcf(raw, name_of_vcf = names_of_vcfs[i]))

    out <- list()
    if (needs_sbs) {
      ann <- annotate_sbs_or_dbs_vcf(sp$SBS, ref_genome = ref_genome,
                                     trans_ranges = trans_ranges,
                                     name_of_vcf = names_of_vcfs[i])
      for (t in intersect(types, c("SBS96", "SBS192", "SBS1536"))) {
        out[[t]] <- vcf_to_sbs_catalog(ann, type = t,
                                       ref_genome = ref_genome,
                                       region = region,
                                       sample_name = names_of_vcfs[i])
      }
    }
    if (needs_dbs || needs_id) {
      stop(
        "DBS and ID catalog types are not yet wired through ",
        "vcfs_to_catalogs() in this port."
      )
    }
    out
  })

  # Column-bind each type across samples
  out <- list()
  for (t in types) {
    cols <- lapply(per_sample, function(s) s[[t]])
    out[[t]] <- cbind_catalogs(cols)
  }
  out
}
