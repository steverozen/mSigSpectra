#' Read transcript ranges from a GENCODE-derived CSV
#'
#' Reads a CSV with columns `chrom`, `start`, `end`, `strand`,
#' `Ensembl.gene.ID`, `gene.symbol` (1-based coordinates), orders the
#' chromosome factor canonically (1..22/19, X, Y), and returns a keyed
#' `data.table`.
#'
#' @param file Path to the transcript-range CSV file.
#'
#' @return A `data.table` keyed on `chrom`, `start`, `end`.
#'
#' @keywords internal
read_transcript_ranges <- function(file) {
  dt <- data.table::fread(file)
  data.table::setnames(dt, c(
    "chrom", "start", "end", "strand", "Ensembl.gene.ID", "gene.symbol"
  ))

  n_chroms <- length(unique(dt$chrom))
  chr_order <- if (n_chroms == 24L) {
    c(as.character(1:22), "X", "Y")
  } else if (n_chroms == 21L) {
    c(as.character(1:19), "X", "Y")
  } else {
    sort(unique(as.character(dt$chrom)))
  }

  dt$chrom <- factor(dt$chrom, chr_order, ordered = TRUE)
  data.table::setkeyv(dt, c("chrom", "start", "end"))
  dt
}

#' Infer transcript ranges for a reference genome
#'
#' If `trans_ranges` is supplied, it is returned unchanged. Otherwise,
#' returns the shipped transcript-ranges table for the given reference
#' genome (GRCh37 / GRCh38 / GRCm38).
#'
#' @param ref_genome A BSgenome object or a character identifier accepted by
#'   [normalize_genome_arg()].
#' @param trans_ranges Optional user-supplied transcript ranges.
#'
#' @return A `data.table` of transcript ranges, or `NULL` if no shipped
#'   table is available and none was supplied.
#'
#' @export
infer_trans_ranges <- function(ref_genome, trans_ranges = NULL) {
  if (!is.null(trans_ranges)) return(trans_ranges)
  if (is_grch37(ref_genome)) return(mSigSpectra::trans.ranges.GRCh37)
  if (is_grch38(ref_genome)) return(mSigSpectra::trans.ranges.GRCh38)
  if (is_grcm38(ref_genome)) return(mSigSpectra::trans.ranges.GRCm38)
  NULL
}

#' Is this reference genome GRCh37 (1000 Genomes hs37d5)?
#'
#' @param x A BSgenome object or a character identifier.
#'
#' @keywords internal
is_grch37 <- function(x) {
  if (is.null(x)) return(FALSE)
  normalize_genome_arg(x)@pkgname == "BSgenome.Hsapiens.1000genomes.hs37d5"
}

#' Is this reference genome GRCh38 (UCSC hg38)?
#'
#' @param x A BSgenome object or a character identifier.
#'
#' @keywords internal
is_grch38 <- function(x) {
  if (is.null(x)) return(FALSE)
  normalize_genome_arg(x)@pkgname == "BSgenome.Hsapiens.UCSC.hg38"
}

#' Is this reference genome GRCm38 (UCSC mm10)?
#'
#' @param x A BSgenome object or a character identifier.
#'
#' @keywords internal
is_grcm38 <- function(x) {
  if (is.null(x)) return(FALSE)
  normalize_genome_arg(x)@pkgname == "BSgenome.Mmusculus.UCSC.mm10"
}

# Defensive column-name collision handlers: if the input VCF already has
# columns that we are about to add, rename the originals with an `_old`
# suffix and warn.

rename_column_if_present <- function(df, col) {
  if (col %in% colnames(df)) {
    new_name <- paste0(col, "_old")
    colnames(df)[colnames(df) == col] <- new_name
    warning(
      'A column named "', col, '" in the VCF was renamed to "',
      new_name, '" to avoid conflict with a newly added column named "',
      col, '".'
    )
  }
  df
}

#' Annotate a VCF data frame with transcript strand information
#'
#' For each variant, finds overlapping transcripts in `trans_ranges` via
#' [GenomicRanges::findOverlaps()] with `type = "within"`, and appends
#' columns `trans.start.pos`, `trans.end.pos`, `trans.strand`,
#' `trans.Ensembl.gene.ID`, `trans.gene.symbol`, plus `bothstrand` (TRUE if
#' the variant falls on transcripts from both strands) and `count` (number
#' of overlapping transcripts).
#'
#' If `trans_ranges` is `NULL` and no shipped table is available for the
#' given `ref_genome`, returns `df` as a `data.table` unchanged (no strand
#' columns added). If `df` has zero rows, it is returned unchanged.
#'
#' @param df A VCF as a data frame / data.table with columns `CHROM`,
#'   `POS`, `ALT`.
#' @param ref_genome A BSgenome object or a character identifier accepted by
#'   [normalize_genome_arg()].
#' @param trans_ranges Optional `data.table` of transcript ranges with
#'   columns `chrom`, `start`, `end`, `strand`, `Ensembl.gene.ID`,
#'   `gene.symbol`. If `NULL`, the shipped table for `ref_genome` is used.
#' @param name_of_vcf Optional VCF name used in warning/error messages.
#'
#' @return A `data.table` with the annotation columns added (left-join
#'   semantics: variants outside any transcript have `NA` values).
#'
#' @export
add_transcript_strand <- function(df, ref_genome, trans_ranges = NULL,
                                  name_of_vcf = NULL) {
  if (nrow(df) == 0L) return(df)

  df$CHROM <- as.character(df$CHROM)

  trans_ranges <- infer_trans_ranges(ref_genome, trans_ranges)
  if (is.null(trans_ranges)) return(data.table::as.data.table(df))

  for (col in c("strand", "start", "end")) {
    df <- rename_column_if_present(df, col)
  }

  ref_genome <- normalize_genome_arg(ref_genome)

  new_chr_names <- check_and_fix_chrom_names_for_trans_ranges(
    trans.ranges = trans_ranges,
    vcf.df = df,
    ref.genome = ref_genome,
    name.of.VCF = name_of_vcf
  )
  trans_ranges <- data.table::copy(trans_ranges)
  trans_ranges$chrom <- new_chr_names

  query_gr <- GenomicRanges::GRanges(
    seqnames = df$CHROM,
    ranges = IRanges::IRanges(start = df$POS, end = df$POS)
  )
  subject_gr <- GenomicRanges::GRanges(
    seqnames = trans_ranges$chrom,
    ranges = IRanges::IRanges(start = trans_ranges$start,
                              end = trans_ranges$end)
  )
  hits <- GenomicRanges::findOverlaps(query_gr, subject_gr, type = "within")
  qi <- S4Vectors::queryHits(hits)
  si <- S4Vectors::subjectHits(hits)

  matched_dt <- data.table::as.data.table(df[qi, , drop = FALSE])
  matched_dt[, `:=`(
    .orig_row       = qi,
    POS2            = df$POS[qi],
    chrom.y         = as.character(trans_ranges$chrom[si]),
    start           = trans_ranges$start[si],
    end             = trans_ranges$end[si],
    strand          = trans_ranges$strand[si],
    Ensembl.gene.ID = trans_ranges$Ensembl.gene.ID[si],
    gene.symbol     = trans_ranges$gene.symbol[si]
  )]

  unmatched_idx <- setdiff(seq_len(nrow(df)), qi)
  if (length(unmatched_idx) > 0L) {
    unmatched_dt <- data.table::as.data.table(df[unmatched_idx, , drop = FALSE])
    unmatched_dt[, `:=`(.orig_row = unmatched_idx, POS2 = df$POS[unmatched_idx])]
    dt <- data.table::rbindlist(list(matched_dt, unmatched_dt), fill = TRUE)
  } else {
    dt <- matched_dt
  }
  data.table::setorder(dt, .orig_row)
  dt[, .orig_row := NULL]

  dt[, bothstrand := ("+" %in% strand) && ("-" %in% strand),
     by = .(CHROM, ALT, POS)]
  dt[, count := .N, by = .(CHROM, ALT, POS)]

  dt[strand == "-", c("end", "start") := .(start, end)]

  df_colnames <- colnames(df)
  trans_colnames <- setdiff(
    c("start", "end", "strand", "Ensembl.gene.ID", "gene.symbol"),
    df_colnames
  )
  data.table::setcolorder(
    dt,
    neworder = c(df_colnames, trans_colnames, "POS2", "chrom.y",
                 "bothstrand", "count")
  )

  data.table::setnames(
    dt,
    old = c("start", "end", "strand", "Ensembl.gene.ID", "gene.symbol"),
    new = c("trans.start.pos", "trans.end.pos", "trans.strand",
            "trans.Ensembl.gene.ID", "trans.gene.symbol")
  )

  dt[, chrom.y := NULL]
  dt
}
