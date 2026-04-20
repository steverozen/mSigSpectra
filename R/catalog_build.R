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

  stop(
    "vcf_to_catalog(type = '", type, "') is not yet implemented in this port. ",
    "Currently only SBS96 / SBS192 / SBS1536 are supported."
  )
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

#' Convenience wrapper: read, annotate, and catalog a set of VCF files
#'
#' Thin composition of [read_vcf()] + [split_vcf()] + [annotate_vcf()] +
#' [vcf_to_catalog()] for a list of files.
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
      ann <- annotate_vcf(sp$SBS, ref_genome = ref_genome,
                          variant_type = "SBS",
                          trans_ranges = trans_ranges,
                          name_of_vcf = names_of_vcfs[i])
      for (t in intersect(types, c("SBS96", "SBS192", "SBS1536"))) {
        out[[t]] <- vcf_to_catalog(ann, type = t, ref_genome = ref_genome,
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
