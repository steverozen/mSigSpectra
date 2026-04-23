#' Build a mutational-spectrum catalog from an annotated VCF
#'
#' Dispatches to a type-specific internal builder based on `type`.
#' Supports SBS96 / SBS192 / SBS1536 today; DBS and ID catalog types are
#' stubs pending the DBS- and indel-specific ports.
#'
#' @param annotated_vcf A VCF annotated by [annotate_vcf()]. For SBS
#'   types, must contain a `seq.<N>bases` column and (for SBS192) a
#'   `trans.strand` / `bothstrand` column.
#' @param type Catalog type identifier (one of `"SBS96"`, `"SBS192"`,
#'   `"SBS1536"`, `"DBS78"`, `"DBS136"`, `"DBS144"`, `"ID83"`, `"ID89"`,
#'   `"ID166"`, `"ID476"`).
#' @param ref_genome Optional BSgenome object or alias; recorded on the
#'   output catalog.
#' @param region One of `"genome"`, `"exome"`, `"transcript"`,
#'   `"unknown"`.
#' @param sample_name Column name for the single-sample catalog matrix
#'   (default `"count"`).
#'
#' @return A single-column numeric matrix with catalog attributes (see
#'   [as_catalog()]).
#'
#' @export
vcf_to_catalog <- function(annotated_vcf,
                           type = c("SBS96", "SBS192", "SBS1536",
                                    "DBS78", "DBS136", "DBS144",
                                    "ID83", "ID89", "ID166", "ID476"),
                           ref_genome = NULL,
                           region = "unknown",
                           sample_name = "count") {
  type <- match.arg(type)
  stop_if_region_illegal(region)

  if (type %in% c("SBS96", "SBS192", "SBS1536")) {
    return(vcf_to_sbs_catalog(annotated_vcf, type, ref_genome, region,
                              sample_name))
  }
  if (type %in% c("DBS78", "DBS136", "DBS144")) {
    return(vcf_to_dbs_catalog(annotated_vcf, type, ref_genome, region,
                              sample_name))
  }
  if (type %in% c("ID83", "ID89", "ID476")) {
    return(vcf_to_id_catalog(annotated_vcf, type, ref_genome, region,
                             sample_name))
  }

  stop(
    "vcf_to_catalog(type = '", type, "') is not yet implemented. ",
    "Supported types: SBS96/192/1536, DBS78/136/144, ID83/89/476."
  )
}

# -----------------------------------------------------------------------------
# DBS catalog builders
# -----------------------------------------------------------------------------

vcf_to_dbs_catalog <- function(vcf, type, ref_genome, region, sample_name) {
  if (nrow(vcf) == 0L) {
    return(empty_dbs_catalog(type, ref_genome, region, sample_name))
  }

  seq_col <- find_seq_context_col(vcf)
  if (is.null(seq_col)) {
    stop("Input VCF has no seq.<N>bases column; run annotate_vcf() first.")
  }
  half_width <- (nchar(vcf[[seq_col]][1]) - 1L) %/% 2L
  center_pos <- half_width + 1L
  quad_start <- center_pos - 1L
  quad_end   <- center_pos + 2L

  stopifnot(all(nchar(vcf$ALT) == 2L))
  stopifnot(all(nchar(vcf$REF) == 2L))

  # Drop rows with N in the tetranucleotide context
  has_n <- grep("N", substr(vcf[[seq_col]], quad_start, quad_end))
  if (length(has_n) > 0L) {
    warning(
      "In sample ", sample_name, ", ", length(has_n),
      " DBS variants have tetranucleotide contexts containing N and were discarded."
    )
    vcf <- vcf[-has_n, ]
  }
  if (nrow(vcf) == 0L) {
    return(empty_dbs_catalog(type, ref_genome, region, sample_name))
  }

  # Deduplicate by (CHROM, ALT, POS) -- one DBS can appear on multiple overlapping
  # transcripts after annotation.
  dedup <- unique(vcf, by = c("CHROM", "ALT", "POS"))

  mat <- NULL
  if (type == "DBS78") {
    canon <- canonicalize_dbs(dedup$REF, dedup$ALT)
    mat <- tabulate_to_catalog_matrix(
      canon, mSigSpectra::catalog.row.order$DBS78, sample_name
    )
  } else if (type == "DBS136") {
    quad <- canonicalize_quad(substr(dedup[[seq_col]], quad_start, quad_end))
    mat <- tabulate_to_catalog_matrix(
      quad, mSigSpectra::catalog.row.order$DBS136, sample_name
    )
  } else if (type == "DBS144") {
    if (!all(c("trans.strand", "bothstrand") %in% colnames(dedup))) {
      stop(
        "vcf_to_catalog(type = 'DBS144') requires trans.strand + bothstrand; ",
        "run annotate_vcf() with a valid trans_ranges first."
      )
    }
    keep <- !is.na(dedup$trans.strand) & (dedup$bothstrand == FALSE)
    sub <- dedup[keep, ]
    if (nrow(sub) == 0L) {
      return(empty_dbs_catalog(type, ref_genome, region, sample_name))
    }
    labels <- paste0(sub$REF, sub$ALT)
    rev_labels <- revc_dbs144(labels)
    labels_final <- ifelse(sub$trans.strand == "-", rev_labels, labels)
    mat <- tabulate_to_catalog_matrix(
      labels_final, mSigSpectra::catalog.row.order$DBS144, sample_name
    )
  }

  as_catalog(mat, type = type, ref_genome = ref_genome, region = region)
}

# Reverse-complement DBS78 REF+ALT pairs whose orientation is non-canonical.
canonicalize_dbs <- function(ref_vec, alt_vec) {
  dbs <- paste0(ref_vec, alt_vec)
  idx <- which(!(dbs %in% mSigSpectra::catalog.row.order$DBS78))
  if (length(idx) == 0L) return(dbs)
  dbs[idx] <- paste0(revc(ref_vec[idx]), revc(alt_vec[idx]))
  stopifnot(all(dbs %in% mSigSpectra::catalog.row.order$DBS78))
  dbs
}

# Reverse-complement DBS136 tetranucleotides whose orientation is non-canonical.
canonicalize_quad <- function(quad) {
  idx <- which(!(quad %in% mSigSpectra::catalog.row.order$DBS136))
  if (length(idx) == 0L) return(quad)
  quad[idx] <- revc(quad[idx])
  stopifnot(all(quad %in% mSigSpectra::catalog.row.order$DBS136))
  quad
}

empty_dbs_catalog <- function(type, ref_genome, region, sample_name) {
  rns <- mSigSpectra::catalog.row.order[[type]]
  m <- matrix(0, nrow = length(rns), ncol = 1L,
              dimnames = list(rns, sample_name))
  as_catalog(m, type = type, ref_genome = ref_genome, region = region)
}

# -----------------------------------------------------------------------------
# ID catalog builder
# -----------------------------------------------------------------------------

# Turn an indel-annotated VCF (with COSMIC_83 / Koh_89 / Koh_476 columns as
# produced by annotate_vcf(variant_type = "ID") + categorize_indels_in_vcf)
# into a count matrix for the requested ID classification scheme.
vcf_to_id_catalog <- function(vcf, type, ref_genome, region, sample_name) {
  if (nrow(vcf) == 0L) {
    rns <- mSigSpectra::catalog.row.order[[if (type == "ID83") "ID" else type]]
    m <- matrix(0, nrow = length(rns), ncol = 1L,
                dimnames = list(rns, sample_name))
    return(as_catalog(m, type = type, ref_genome = ref_genome, region = region))
  }

  col <- switch(type, ID83 = "COSMIC_83", ID89 = "Koh_89", ID476 = "Koh_476")
  if (!col %in% colnames(vcf)) {
    stop(
      "vcf_to_catalog(type = '", type, "') requires column '", col,
      "' on the input; run annotate_vcf(variant_type = 'ID') first."
    )
  }

  df <- as.data.frame(vcf)
  res <- switch(
    type,
    ID83  = annot_vcf_to_83_catalog(df,  sample_id = sample_name),
    ID89  = annot_vcf_to_89_catalog(df,  sample_id = sample_name),
    ID476 = annot_vcf_to_476_catalog(df, sample_id = sample_name)
  )
  m <- as.matrix(res)
  storage.mode(m) <- "numeric"
  as_catalog(m, type = type, ref_genome = ref_genome, region = region)
}

# Build SBS96 / SBS192 / SBS1536 from an annotated SBS VCF. Returns a
# single catalog matrix of the requested `type`; intermediate matrices
# for the other SBS resolutions are still computed (cheap) but not
# returned -- this keeps the public API focused on "one call, one
# catalog type".

vcf_to_sbs_catalog <- function(vcf, type, ref_genome, region, sample_name) {
  seq_col <- find_seq_context_col(vcf)
  if (is.null(seq_col) && nrow(vcf) > 0L) {
    stop(
      "Input VCF has no seq.<N>bases column; run annotate_vcf() first ",
      "with variant_type = 'SBS'."
    )
  }

  # Short-circuit for empty input.
  if (nrow(vcf) == 0L) {
    return(empty_sbs_catalog(type, ref_genome, region, sample_name))
  }

  stopifnot(all(nchar(vcf$ALT) == 1L))
  stopifnot(all(nchar(vcf$REF) == 1L))
  stopifnot(all(vcf$ALT != vcf$REF))

  # Drop rows where the REF allele doesn't match the reference genome
  # pentanucleotide (VCF / genome mismatch), or where the pentanucleotide
  # context contains N.
  half_width <- (nchar(vcf[[seq_col]][1]) - 1L) %/% 2L
  center_pos <- half_width + 1L
  penta_start <- center_pos - 2L
  penta_end   <- center_pos + 2L

  mismatches <- which(
    vcf$REF != substr(vcf[[seq_col]], center_pos, center_pos)
  )
  if (length(mismatches) > 0L) {
    message(
      "In sample ", sample_name, ", ", length(mismatches), " of ",
      nrow(vcf), " SBS variants have REF != reference-genome base and ",
      "have been discarded."
    )
    vcf <- vcf[-mismatches, ]
  }
  has_n <- grep("N", substr(vcf[[seq_col]], penta_start, penta_end))
  if (length(has_n) > 0L) {
    warning(
      "In sample ", sample_name, ", ", length(has_n), " SBS variants have ",
      "pentanucleotide contexts containing N and have been discarded."
    )
    vcf <- vcf[-has_n, ]
  }

  if (nrow(vcf) == 0L) {
    return(empty_sbs_catalog(type, ref_genome, region, sample_name))
  }

  # Build pentanucleotide + ALT string, canonicalize to pyrimidine form.
  context <- substr(vcf[[seq_col]], penta_start, penta_end)
  vcf$mutation <- paste0(context, vcf$ALT)
  vcf$pyr.mut  <- pyr_penta(vcf$mutation)

  # One SBS can appear on multiple transcripts. Deduplicate on
  # (CHROM, ALT, POS) before counting SBS96 / SBS1536.
  dedup <- unique(vcf, by = c("CHROM", "ALT", "POS"))

  # SBS1536 counts
  mat1536 <- tabulate_to_catalog_matrix(
    dedup$pyr.mut,
    row_order = mSigSpectra::catalog.row.order$SBS1536,
    sample_name = sample_name
  )

  # SBS96 counts derived from SBS1536 (collapse 5mer -> 3mer + ALT)
  # Row label for SBS96: central trinucleotide + ALT, i.e.
  # substr(pyr.mut, 2, 4) + substr(pyr.mut, 6, 6).
  sbs96_rows <- paste0(
    substr(rownames(mat1536), 2L, 4L),
    substr(rownames(mat1536), 6L, 6L)
  )
  mat96_rawdt <- data.table::data.table(rn = sbs96_rows, count = mat1536[, 1L])
  mat96_dt <- mat96_rawdt[, .(count = sum(count)), by = rn]
  mat96 <- matrix(
    mat96_dt$count[match(mSigSpectra::catalog.row.order$SBS96, mat96_dt$rn)],
    ncol = 1L,
    dimnames = list(mSigSpectra::catalog.row.order$SBS96, sample_name)
  )
  mat96[is.na(mat96)] <- 0

  # SBS192 (stranded) — requires trans.strand / bothstrand columns.
  mat192 <- NULL
  if (type == "SBS192") {
    if (!all(c("trans.strand", "bothstrand") %in% colnames(vcf))) {
      stop(
        "vcf_to_catalog(type = 'SBS192') requires trans.strand / bothstrand ",
        "columns on the input; run annotate_vcf() with a non-NULL ",
        "trans_ranges or a supported ref_genome first."
      )
    }
    mat192 <- build_sbs192_matrix(vcf, sample_name)
  }

  out <- switch(
    type,
    SBS96   = mat96,
    SBS1536 = mat1536,
    SBS192  = mat192
  )

  as_catalog(out, type = type, ref_genome = ref_genome, region = region)
}

find_seq_context_col <- function(vcf) {
  cols <- grep("^seq\\.\\d+bases$", colnames(vcf), value = TRUE)
  if (length(cols) == 0L) return(NULL)
  cols[1L]
}

tabulate_to_catalog_matrix <- function(values, row_order, sample_name) {
  tab <- table(factor(values, levels = row_order))
  m <- matrix(as.numeric(tab), ncol = 1L,
              dimnames = list(names(tab), sample_name))
  m[row_order, , drop = FALSE]
}

empty_sbs_catalog <- function(type, ref_genome, region, sample_name) {
  key <- if (type == "ID83") "ID" else type
  rns <- mSigSpectra::catalog.row.order[[key]]
  m <- matrix(0, nrow = length(rns), ncol = 1L,
              dimnames = list(rns, sample_name))
  as_catalog(m, type = type, ref_genome = ref_genome, region = region)
}

build_sbs192_matrix <- function(vcf, sample_name) {
  # Only count rows that have a non-NA strand and are not on bothstrand
  keep <- !is.na(vcf$trans.strand) & (vcf$bothstrand == FALSE)
  sub <- vcf[keep, ]

  # Deduplicate on (CHROM, ALT, POS) within transcribed subset
  sub <- unique(sub, by = c("CHROM", "ALT", "POS"))

  if (nrow(sub) == 0L) {
    rns <- mSigSpectra::catalog.row.order$SBS192
    return(matrix(0, nrow = length(rns), ncol = 1L,
                  dimnames = list(rns, sample_name)))
  }

  # SBS192 row label = trinucleotide + ALT, oriented by transcript strand:
  #   '+' strand: use pyr-canonical form directly.
  #   '-' strand: reverse-complement to put onto the pyr-canonical strand.
  # We already computed pyr.mut on the unstranded strand; for strand-aware
  # counting we need the mutation on the transcript strand itself.
  trin <- substr(sub$mutation, 2L, 4L)
  alt_base <- substr(sub$mutation, 6L, 6L)
  labels <- paste0(trin, alt_base)
  labels_rev <- revc_sbs96(labels)
  labels_final <- ifelse(sub$trans.strand == "-", labels_rev, labels)

  dt <- data.table::data.table(rn = labels_final)
  dt_counts <- dt[, .(count = .N), by = rn]

  rns <- mSigSpectra::catalog.row.order$SBS192
  m <- matrix(
    dt_counts$count[match(rns, dt_counts$rn)],
    ncol = 1L,
    dimnames = list(rns, sample_name)
  )
  m[is.na(m)] <- 0
  m
}

