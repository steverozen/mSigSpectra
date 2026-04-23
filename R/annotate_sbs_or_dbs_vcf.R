#' Annotate an SBS or DBS VCF with flanking sequence context and transcript
#' strand
#'
#' Adds flanking sequence context via [add_seq_context()] and, when a
#' transcript-ranges table is available, transcript strand via
#' [add_transcript_strand()]. The same pipeline is appropriate for SBS and
#' DBS; DBS-specific validation (e.g. no `N` in the tetranucleotide
#' context) happens downstream in [vcf_to_dbs_catalog()].
#'
#' @param vcf A VCF as a data.frame / data.table with `CHROM`, `POS`,
#'   `REF`, `ALT` columns. For SBS, `REF` and `ALT` are single bases; for
#'   DBS, both are two bases.
#' @param ref_genome A BSgenome object or a character alias accepted by
#'   [normalize_genome_arg()].
#' @param trans_ranges Optional transcript-ranges `data.table`. If `NULL`,
#'   the shipped table is used when available for `ref_genome`.
#' @param seq_context_width Width (per side) of the flanking-sequence
#'   window (default 10 → 21-base window).
#' @param name_of_vcf Optional VCF name used in warnings.
#'
#' @return A list with:
#'   * `annotated.vcf`: the input VCF with new columns `seq.<N>bases` and
#'     (when transcript ranges are available) `trans.strand` /
#'     `bothstrand`.
#'   * `discarded.variants`: currently always `NULL` for this path
#'     (included for symmetry with [annotate_id_vcf()]).
#'
#' @export
annotate_sbs_or_dbs_vcf <- function(vcf,
                                    ref_genome,
                                    trans_ranges = NULL,
                                    seq_context_width = 10L,
                                    name_of_vcf = NULL) {
  vcf <- data.table::as.data.table(vcf)
  if (nrow(vcf) == 0L) {
    return(list(annotated.vcf = vcf, discarded.variants = NULL))
  }

  vcf <- add_seq_context(
    vcf, ref_genome = ref_genome,
    seq_context_width = seq_context_width, name_of_vcf = name_of_vcf
  )

  resolved_ranges <- infer_trans_ranges(ref_genome, trans_ranges)
  if (!is.null(resolved_ranges)) {
    vcf <- add_transcript_strand(
      vcf, ref_genome = ref_genome,
      trans_ranges = resolved_ranges, name_of_vcf = name_of_vcf
    )
  }

  list(
    annotated.vcf = data.table::as.data.table(vcf),
    discarded.variants = NULL
  )
}
