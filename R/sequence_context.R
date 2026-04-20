#' Add flanking sequence context to a VCF data frame
#'
#' Extracts `seq_context_width` bases upstream and downstream of each variant
#' position from the reference genome and attaches them as a new column
#' named `seq.<N>bases` where `N = 2 * seq_context_width + 1`.
#'
#' @param df A VCF as a data frame / data.table with columns `CHROM` and
#'   `POS`.
#' @param ref_genome A BSgenome object or a character identifier accepted by
#'   [normalize_genome_arg()].
#' @param seq_context_width Number of flanking bases on each side (default 10,
#'   producing a 21-base window).
#' @param name_of_vcf Optional VCF name used in warning/error messages.
#'
#' @return `df` with an added character column `seq.<N>bases`.
#'
#' @export
add_seq_context <- function(df, ref_genome, seq_context_width = 10,
                            name_of_vcf = NULL) {
  if (nrow(df) == 0L) return(df)

  ref_genome <- normalize_genome_arg(ref_genome)

  chr_names <- check_and_fix_chrom_names(
    vcf.df = df,
    ref.genome = ref_genome,
    name.of.VCF = name_of_vcf
  )

  ranges <- GenomicRanges::GRanges(
    chr_names,
    IRanges::IRanges(
      start = df$POS - seq_context_width,
      end   = df$POS + seq_context_width
    )
  )

  extracted <- BSgenome::getSeq(ref_genome, ranges, as.character = TRUE)
  col <- paste0("seq.", 2L * seq_context_width + 1L, "bases")
  df[[col]] <- extracted
  df
}
