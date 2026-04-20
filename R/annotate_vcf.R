#' Annotate a VCF with flanking sequence context and transcript strand
#'
#' `annotate_vcf()` is the high-level annotator used by
#' [vcf_to_catalog()]. It dispatches on `variant_type`:
#'
#' * `"SBS"`: add flanking sequence context ([add_seq_context()]) and
#'   transcript strand ([add_transcript_strand()] if `trans_ranges` is
#'   supplied or shipped for `ref_genome`).
#' * `"DBS"`: same as SBS; downstream [vcf_to_catalog()] will require the
#'   tetranucleotide context to be free of `N`.
#' * `"ID"`: placeholder — indel annotation is handled by a separate
#'   module ([annotate_id_vcf()]) that additionally runs indel
#'   justification and repeat-context detection. Not yet implemented in
#'   this port.
#' * `"auto"`: classify each row by REF/ALT length and annotate each
#'   subset appropriately. Only the SBS / DBS paths are wired up for now.
#'
#' @param vcf A VCF as a data.frame / data.table with `CHROM`, `POS`,
#'   `REF`, `ALT` columns.
#' @param ref_genome A BSgenome object or a character alias accepted by
#'   [normalize_genome_arg()].
#' @param variant_type One of `"auto"`, `"SBS"`, `"DBS"`, `"ID"`. For
#'   `"auto"`, rows are routed by REF/ALT length.
#' @param trans_ranges Optional transcript-ranges `data.table`. If `NULL`,
#'   the shipped table is used when available.
#' @param seq_context_width Width (per side) of the flanking-sequence
#'   window; default 10 (produces a 21-base window).
#' @param name_of_vcf Optional VCF name used in warnings.
#'
#' @return A `data.table` with the annotation columns added.
#'
#' @export
annotate_vcf <- function(vcf, ref_genome,
                         variant_type = c("auto", "SBS", "DBS", "ID"),
                         trans_ranges = NULL,
                         seq_context_width = 10L,
                         name_of_vcf = NULL) {
  variant_type <- match.arg(variant_type)
  vcf <- data.table::as.data.table(vcf)
  if (nrow(vcf) == 0L) return(vcf)

  if (variant_type == "auto") {
    ref_n <- nchar(vcf$REF)
    alt_n <- nchar(vcf$ALT)
    inferred <- if (all(ref_n == 1L & alt_n == 1L)) {
      "SBS"
    } else if (all(ref_n == 2L & alt_n == 2L)) {
      "DBS"
    } else if (all(ref_n != alt_n)) {
      "ID"
    } else {
      stop(
        "variant_type='auto' requires a homogeneous VCF; ",
        "this VCF has a mix of SBS/DBS/ID rows. ",
        "Use split_vcf() to partition first, or pass variant_type explicitly."
      )
    }
    variant_type <- inferred
  }

  switch(
    variant_type,
    SBS = annotate_sbs_vcf(vcf, ref_genome, trans_ranges,
                           seq_context_width, name_of_vcf),
    DBS = annotate_dbs_vcf(vcf, ref_genome, trans_ranges,
                           seq_context_width, name_of_vcf),
    ID  = stop(
      "annotate_vcf(variant_type = 'ID') is not yet implemented in this port. ",
      "Indel annotation requires the indel_justify / indel_categorize modules."
    )
  )
}

annotate_sbs_vcf <- function(vcf, ref_genome, trans_ranges,
                             seq_context_width, name_of_vcf) {
  vcf <- add_seq_context(
    vcf, ref_genome = ref_genome,
    seq_context_width = seq_context_width, name_of_vcf = name_of_vcf
  )

  trans_ranges <- infer_trans_ranges(ref_genome, trans_ranges)
  if (!is.null(trans_ranges)) {
    vcf <- add_transcript_strand(
      vcf, ref_genome = ref_genome,
      trans_ranges = trans_ranges, name_of_vcf = name_of_vcf
    )
  }
  data.table::as.data.table(vcf)
}

annotate_dbs_vcf <- function(vcf, ref_genome, trans_ranges,
                             seq_context_width, name_of_vcf) {
  # Same annotation pipeline as SBS. DBS-specific validation (no N in
  # tetranucleotide context) is done by vcf_to_catalog() on the way to the
  # DBS catalog.
  annotate_sbs_vcf(vcf, ref_genome, trans_ranges,
                   seq_context_width, name_of_vcf)
}
