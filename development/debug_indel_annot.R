# Debug / exploration scaffolding for indel annotation and categorization.
# Run interactively after devtools::load_all().
#
# Available test VCFs (inst/extdata == tests/testthat/testdata, files are identical):
#   Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf   408 indel variants
#   Mutect-GRCh37/Mutect.GRCh37.s1.vcf            (SBS; not indels)
#   Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf (SBS; not indels)

# ---- helpers ----------------------------------------------------------------

#' Read the bundled Strelka indel VCF and return an in-memory data.frame.
#'
#' @param n Maximum number of variants to return (NULL = all).
#' @param filter Passed to read_vcf(); "PASS" to keep only PASS variants.
#' @param from_inst If TRUE use inst/extdata; if FALSE use tests/testthat/testdata.
load_strelka_id_vcf <- function(n = 20L, filter = "PASS", from_inst = TRUE) {
  base <- if (from_inst) {
    system.file(
      "extdata",
      "Strelka-ID-GRCh37",
      "Strelka.ID.GRCh37.s1.vcf",
      package = "mSigSpectra"
    )
  } else {
    testthat::test_path(
      "testdata",
      "Strelka-ID-GRCh37",
      "Strelka.ID.GRCh37.s1.vcf"
    )
  }
  vcf <- read_vcf(base, filter = filter)
  if (!is.null(n)) {
    vcf <- vcf[seq_len(min(n, nrow(vcf))), ]
  }
  vcf
}

#' Justify and fully annotate (seq context + categories) a small indel VCF.
#'
#' @param vcf data.frame from load_strelka_id_vcf() or similar.
#' @param ref_genome Genome alias, e.g. "GRCh37".
load_annotated_id_vcf <- function(vcf, ref_genome = "GRCh37") {
  annotate_id_vcf(vcf, ref_genome = ref_genome, explain_indels = 0)
}

# ---- quick end-to-end run ---------------------------------------------------

if (TRUE) {
  devtools::load_all()

  vcf_raw <- load_strelka_id_vcf(n = 20L)
  result <- load_annotated_id_vcf(vcf_raw)
  ann <- result$annotated.vcf
  discarded <- result$discarded.variants

  # Inspect the columns added by categorize_indels_in_vcf():
  cat_cols <- c(
    "ins_or_del",
    "pre",
    "ins_or_del_seq",
    "post",
    "L",
    "U_seq",
    "U",
    "U_seq_count_in_indel_seq",
    "indel_str_count_in_ref",
    "R",
    "R_outside_ins_or_del_seq",
    "mh",
    "short_visual",
    "long_visual",
    # from koh_extra (seg_simple):
    "unit",
    "unit_length",
    "internal_rep",
    "internal_reps",
    "spacer",
    "spacer_length",
    "prime3_rep",
    "prime3_reps",
    "original_reps",
    # classification strings:
    "COSMIC_83",
    "Koh_89",
    "Koh_476"
  )
  print(ann[, intersect(cat_cols, colnames(ann))])

  # Drive categorize_1_justified_indel() on a single row for debugging:
  row <- ann[1, ]
  cw <- row$seq.context.width
  res <- categorize_1_justified_indel(
    context = row$seq.context,
    ins_or_del = if (nchar(row$ALT) < nchar(row$REF)) "d" else "i",
    ins_or_del_seq = {
      iod <- if (nchar(row$ALT) < nchar(row$REF)) "d" else "i"
      if (iod == "d") {
        substr(row$REF, 2L, nchar(row$REF))
      } else {
        substr(row$ALT, 2L, nchar(row$ALT))
      }
    },
    pos = cw + 2L
  )
  str(res)

  # Run on all 408 variants (need BSgenome.Hsapiens.1000genomes.hs37d5):
  # vcf_all  <- load_strelka_id_vcf(n = NULL)
  # result_all <- load_annotated_id_vcf(vcf_all)
}
