#' Annotate an in-memory ID (indel) VCF with sequence context, transcript
#' strand, and COSMIC / Koh indel categories
#'
#' Port of ICAMS's `AnnotateIDVCF()`. Justifies each indel against the
#' reference genome, optionally annotates transcript strand, and categorizes
#' each justified indel into the COSMIC 83, Koh 89, and Koh 476 schemes.
#'
#' @param vcf An in-memory ID VCF as a `data.frame` / `data.table`. Expects
#'   a "context base" to the left of the indel, e.g. `REF = "ACG"`,
#'   `ALT = "A"` (deletion of `CG`) or `REF = "A"`, `ALT = "ACC"` (insertion
#'   of `CC`).
#' @param ref_genome A BSgenome object or a character alias (e.g.
#'   `"GRCh37"`, `"GRCh38"`, `"mm10"`).
#' @param trans_ranges Optional transcript-ranges `data.table`. If `NULL`
#'   and `add_transcript_ranges = TRUE`, the shipped table is used when
#'   available for `ref_genome`.
#' @param name_of_vcf Optional name for the VCF, used only in warning /
#'   error messages.
#' @param suppress_discarded_variants_warnings If `TRUE`, do not warn when
#'   variants that cannot be processed are discarded.
#' @param explain_indels If `0`, do not explain. If `1`, emit messages for
#'   indels whose position was left-shifted by justification. If `2`, emit
#'   messages for all indels.
#' @param context_width_multiplier Multiplier used to guess how much
#'   flanking sequence (per side) is needed to categorize each indel.
#' @param add_transcript_ranges If `TRUE`, annotate transcript strand when a
#'   transcript-ranges table is available.
#'
#' @return A list with:
#'   * `annotated.vcf`: the input VCF with new columns `seq.context`,
#'     `seq.context.width`, `pos_shift`, `COSMIC_83`, `Koh_89`, `Koh_476`,
#'     and (when transcript ranges are available) `trans.strand` /
#'     `bothstrand`.
#'   * `discarded.variants`: rows that could not be justified, or `NULL` if
#'     none were discarded.
#'
#' @export
annotate_id_vcf <- function(vcf,
                            ref_genome,
                            trans_ranges = NULL,
                            name_of_vcf = NULL,
                            suppress_discarded_variants_warnings = TRUE,
                            explain_indels = 1,
                            context_width_multiplier = 20L,
                            add_transcript_ranges = TRUE) {
  justified <- justify_id_vcf(
    ID.vcf = vcf,
    ref.genome = ref_genome,
    name.of.VCF = name_of_vcf,
    suppress.discarded.variants.warnings = suppress_discarded_variants_warnings,
    explain_indels = explain_indels,
    context_width_multiplier = context_width_multiplier
  )

  ann <- data.table::as.data.table(justified$annotated.vcf)
  discarded <- justified$discarded.variants

  if (nrow(ann) == 0L) {
    return(list(annotated.vcf = ann, discarded.variants = discarded))
  }

  if (add_transcript_ranges) {
    resolved_ranges <- infer_trans_ranges(ref_genome, trans_ranges)
    if (!is.null(resolved_ranges)) {
      ann <- add_transcript_strand(
        ann,
        ref_genome = ref_genome,
        trans_ranges = resolved_ranges,
        name_of_vcf = name_of_vcf
      )
      ann <- data.table::as.data.table(ann)
    }
  }

  cats <- categorize_indels_in_vcf(ann)
  indel_info_df <- data.table::rbindlist(cats, fill = TRUE)
  ann <- cbind(ann, indel_info_df)

  list(annotated.vcf = ann, discarded.variants = discarded)
}
