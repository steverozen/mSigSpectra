# Row-wise categorizer: for each justified indel in `vcf`, call
# categorize_1_justified_indel() and return a list of per-row result
# lists. Callers flatten via data.table::rbindlist(..., fill = TRUE) to
# add per-indel metadata columns (COSMIC_83, Koh_89, Koh_476, R, L, U,
# ins_or_del, pre, post, mh, ...) to the VCF. Mirrors the shape of
# ICAMS's categorize_indels_in_vcf().
categorize_indels_in_vcf <- function(vcf) {
  n <- nrow(vcf)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    ref <- vcf$REF[i]
    alt <- vcf$ALT[i]
    ctx <- vcf$seq.context[i]
    cw  <- vcf$seq.context.width[i]
    ins_or_del <- if (nchar(alt) < nchar(ref)) "d" else "i"
    ins_or_del_seq <- if (ins_or_del == "d") {
      substr(ref, 2L, nchar(ref))
    } else {
      substr(alt, 2L, nchar(alt))
    }
    pos_in_ctx <- cw + 2L

    out[[i]] <- tryCatch(
      categorize_1_justified_indel(
        context = ctx,
        ins_or_del = ins_or_del,
        ins_or_del_seq = ins_or_del_seq,
        pos = pos_in_ctx
      ),
      error = function(e) list(
        COSMIC_83 = NA_character_,
        Koh_89    = NA_character_,
        Koh_476   = NA_character_
      )
    )
  }
  out
}
