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
#' * `"ID"`: placeholder — indel annotation will require indel
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
    ID  = annotate_id_vcf(vcf, ref_genome, trans_ranges,
                          seq_context_width, name_of_vcf)
  )
}

annotate_id_vcf <- function(vcf, ref_genome, trans_ranges,
                            seq_context_width, name_of_vcf) {
  # 1) Justify each indel (left-shifts canonical position, adds seq.context /
  #    seq.context.width / pos_shift / updated REF, ALT, POS).
  out <- suppressWarnings(justify_id_vcf(
    vcf, ref.genome = ref_genome,
    name.of.VCF = name_of_vcf,
    explain_indels = 0
  ))
  ann <- data.table::as.data.table(out$annotated.vcf)
  if (nrow(ann) == 0L) return(ann)

  # 2) Categorize each justified indel into COSMIC_83 / Koh_89 / Koh_476 strings.
  cats <- categorize_indels_in_vcf(ann)
  ann[, `:=`(
    COSMIC_83 = cats$COSMIC_83,
    Koh_89    = cats$Koh_89,
    Koh_476   = cats$Koh_476
  )]

  # 3) Optional transcript-strand annotation (for stranded ID catalogs, if added
  #    in the future).
  trans_ranges <- infer_trans_ranges(ref_genome, trans_ranges)
  if (!is.null(trans_ranges)) {
    ann <- add_transcript_strand(
      ann, ref_genome = ref_genome,
      trans_ranges = trans_ranges, name_of_vcf = name_of_vcf
    )
  }
  ann
}

# Row-wise categorizer: for each justified indel in `vcf`, call
# categorize_1_justified_indel() and collect the three output strings.
categorize_indels_in_vcf <- function(vcf) {
  n <- nrow(vcf)
  cosmic <- character(n); koh89 <- character(n); koh476 <- character(n)

  for (i in seq_len(n)) {
    ref <- vcf$REF[i]
    alt <- vcf$ALT[i]
    ctx <- vcf$seq.context[i]
    cw  <- vcf$seq.context.width[i]
    # Strip the common first base.
    ins_or_del <- if (nchar(alt) < nchar(ref)) "d" else "i"
    ins_or_del_seq <- if (ins_or_del == "d") {
      substr(ref, 2L, nchar(ref))
    } else {
      substr(alt, 2L, nchar(alt))
    }
    pos_in_ctx <- cw + 2L   # position in seq.context just after the common base

    r <- tryCatch(
      categorize_1_justified_indel(
        context = ctx,
        ins_or_del = ins_or_del,
        ins_or_del_seq = ins_or_del_seq,
        pos = pos_in_ctx
      ),
      error = function(e) list(COSMIC_83 = NA_character_,
                               Koh_89    = NA_character_,
                               Koh_476   = NA_character_)
    )
    cosmic[i] <- r$COSMIC_83 %||% NA_character_
    koh89[i]  <- r$Koh_89    %||% NA_character_
    koh476[i] <- r$Koh_476   %||% NA_character_
  }

  list(COSMIC_83 = cosmic, Koh_89 = koh89, Koh_476 = koh476)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

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
