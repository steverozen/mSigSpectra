#' Transform a catalog between counts and density
#'
#' Converts between counts (raw mutation counts per category) and density
#' (mutations per megabase of context) representations, and between
#' catalog ↔ signature (column-normalized) forms. The transformation is
#' expressed by multiplying each category's count by
#' `target_abundance[source.n.mer] / source_abundance[source.n.mer]`,
#' where `source.n.mer` is the reference context encoded in the row name
#' (first 3 characters for SBS96/192, 5 for SBS1536, 2 for DBS78/144, 4
#' for DBS136).
#'
#' @param catalog An mSigSpectra catalog (see [as_catalog()]) with
#'   attributes `type`, `counts_or_density`, `ref_genome`, `region`,
#'   `abundance`.
#' @param target_ref_genome,target_region,target_counts_or_density Target
#'   attributes. If `NULL`, the catalog's current attribute is used.
#' @param target_abundance Optional target abundance vector. If `NULL`,
#'   inferred from `target_ref_genome` + `target_region`.
#'
#' @return A new catalog (plain matrix with attributes).
#'
#' @export
transform_catalog <- function(catalog,
                              target_ref_genome = NULL,
                              target_region = NULL,
                              target_counts_or_density = NULL,
                              target_abundance = NULL) {
  check_catalog_attributes(catalog, target_counts_or_density)

  s <- list(
    ref_genome        = attr(catalog, "ref_genome",        exact = TRUE),
    counts_or_density = attr(catalog, "counts_or_density", exact = TRUE),
    abundance         = attr(catalog, "abundance",         exact = TRUE),
    region            = attr(catalog, "region",            exact = TRUE),
    type              = attr(catalog, "type",              exact = TRUE)
  )

  if (is.null(target_ref_genome))        target_ref_genome        <- s$ref_genome
  if (is.null(target_counts_or_density)) target_counts_or_density <- s$counts_or_density
  if (is.null(target_region))            target_region            <- s$region
  stop_if_counts_or_density_illegal(target_counts_or_density)
  stop_if_region_illegal(target_region)

  inferred <- infer_abundance(catalog, target_ref_genome,
                              target_region, target_counts_or_density)
  if (is.null(target_abundance)) target_abundance <- inferred

  if (!is.null(s$abundance) && !is.null(target_abundance)) {
    if (!identical(names(s$abundance), names(target_abundance))) {
      stop(
        "The names of target_abundance must be identical to the names ",
        "of the catalog's abundance attribute"
      )
    }
  }

  t <- list(
    ref_genome        = target_ref_genome,
    counts_or_density = target_counts_or_density,
    region            = target_region,
    abundance         = target_abundance,
    type              = s$type
  )

  test <- is_transformation_legal(s, t)
  if (test == "null.op") return(catalog)
  if (test == "sig.only") {
    out <- apply(catalog, 2L, function(x) x / sum(x))
    return(as_catalog(
      out, type = s$type, ref_genome = t$ref_genome,
      region = t$region, abundance = t$abundance,
      counts_or_density = t$counts_or_density
    ))
  }

  factor <- t$abundance / s$abundance
  factor[is.infinite(factor) | is.na(factor)] <- 0
  names(factor) <- names(t$abundance)

  out.catalog <- catalog
  for (src_n_mer in names(s$abundance)) {
    rows <- grep(paste0("^", src_n_mer), rownames(out.catalog))
    if (length(rows) > 0L) {
      out.catalog[rows, ] <- out.catalog[rows, ] * factor[[src_n_mer]]
    }
  }

  if (t$counts_or_density %in% c("counts.signature", "density.signature")) {
    out.catalog <- apply(out.catalog, 2L, function(x) x / sum(x))
  }

  as_catalog(
    out.catalog, type = s$type, ref_genome = t$ref_genome,
    region = t$region, abundance = t$abundance,
    counts_or_density = t$counts_or_density
  )
}

check_catalog_attributes <- function(catalog, target_counts_or_density) {
  ref_genome        <- attr(catalog, "ref_genome",        exact = TRUE)
  counts_or_density <- attr(catalog, "counts_or_density", exact = TRUE)
  abundance         <- attr(catalog, "abundance",         exact = TRUE)
  region            <- attr(catalog, "region",            exact = TRUE)

  if (!is.null(ref_genome) && !inherits(ref_genome, "BSgenome") &&
      identical(ref_genome, "mixed")) {
    stop('Cannot transform a catalog with ref_genome = "mixed"')
  }
  if (is.null(counts_or_density)) {
    stop("Cannot transform a catalog with NULL counts_or_density")
  }
  if (is.null(abundance) &&
      !(counts_or_density == "counts" &&
        identical(target_counts_or_density, "counts.signature"))) {
    stop("Cannot transform a catalog with NULL abundance")
  }
  if (is.null(region)) {
    stop("Cannot transform a catalog with NULL region")
  }
  if (identical(region, "mixed")) {
    stop('Cannot transform a catalog with region = "mixed"')
  }
  invisible(TRUE)
}

is_signature <- function(x) x %in% c("counts.signature", "density.signature")
is_counts    <- function(x) x %in% c("counts", "counts.signature")

is_transformation_legal <- function(s, t) {
  if (is_signature(s$counts_or_density) && !is_signature(t$counts_or_density)) {
    stop("Cannot transform from a signature to a counts or density catalog")
  }
  if (identical(s$type, "COMPOSITE") || identical(t$type, "COMPOSITE")) {
    stop("Transformation of COMPOSITE catalogs not supported")
  }
  if (is_density_type(s$counts_or_density)) {
    transform_from_density(s, t)
  } else {
    transform_from_counts(s, t)
  }
}

is_density_type <- function(x) x %in% c("density", "density.signature")

transform_from_counts <- function(s, t) {
  if (is_signature(s$counts_or_density)) {
    # counts.signature -> counts.signature / density.signature only
    if (identical(t$counts_or_density, "counts.signature")) {
      if (is.null(s$abundance) || is.null(t$abundance)) {
        stop("counts.signature -> counts.signature requires non-NULL abundance")
      }
      if (all(s$abundance == t$abundance)) return("null.op")
      return(TRUE)
    }
    if (identical(t$counts_or_density, "density.signature")) {
      if (is.null(s$abundance)) {
        stop("counts.signature -> density.signature requires non-NULL source abundance")
      }
      return(TRUE)
    }
    stop("Cannot transform counts.signature -> ", t$counts_or_density)
  }

  if (identical(t$counts_or_density, "counts.signature")) {
    if (abundance_is_same(s$abundance, t$abundance)) return("sig.only")
    return(TRUE)
  }
  if (identical(t$counts_or_density, "counts")) {
    if (is.null(s$abundance) || is.null(t$abundance)) {
      stop("counts -> counts requires non-NULL abundance on both sides")
    }
    if (all(s$abundance == t$abundance)) return("null.op")
    return(TRUE)
  }
  if (identical(t$counts_or_density, "density.signature")) {
    if (is.null(s$abundance)) {
      stop("counts -> density.signature requires non-NULL source abundance")
    }
    return(TRUE)
  }
  if (identical(t$counts_or_density, "density")) {
    if (is.null(s$abundance)) {
      stop("counts -> density requires non-NULL source abundance")
    }
    return(TRUE)
  }
  stop("programming error: unexpected counts_or_density ", t$counts_or_density)
}

transform_from_density <- function(s, t) {
  if (is_signature(s$counts_or_density)) {
    if (identical(t$counts_or_density, "density.signature")) return("null.op")
    if (identical(t$counts_or_density, "counts.signature")) {
      if (is.null(t$abundance)) {
        stop("density.signature -> counts.signature requires non-NULL target abundance")
      }
      return(TRUE)
    }
    stop("Cannot transform density.signature -> ", t$counts_or_density)
  }
  if (identical(t$counts_or_density, "density")) return("null.op")
  if (identical(t$counts_or_density, "density.signature")) return("sig.only")
  if (identical(t$counts_or_density, "counts")) {
    if (is.null(t$abundance)) {
      stop("density -> counts requires non-NULL target abundance")
    }
    return(TRUE)
  }
  if (identical(t$counts_or_density, "counts.signature")) {
    if (is.null(t$abundance)) {
      stop("density -> counts.signature requires non-NULL target abundance")
    }
    return(TRUE)
  }
  stop("programming error: unexpected counts_or_density ", t$counts_or_density)
}

abundance_is_same <- function(a1, a2) {
  if (is.null(a1) && is.null(a2))  return(TRUE)
  if (is.null(a1) || is.null(a2))  return(FALSE)
  all(a1 == a2) && identical(names(a1), names(a2))
}

#' Collapse a higher-resolution catalog to a lower-resolution one
#'
#' Supports:
#' * `SBS1536 -> SBS96`: collapse pentanucleotide → trinucleotide contexts.
#' * `SBS192  -> SBS96`: sum across the two transcript strands.
#' * `DBS144  -> DBS78`: sum across the two transcript strands.
#'
#' @param catalog An mSigSpectra catalog.
#' @param to Target catalog type (`"SBS96"`, `"DBS78"`).
#'
#' @return A new catalog with the collapsed rows.
#'
#' @export
collapse_catalog <- function(catalog, to = c("SBS96", "DBS78")) {
  to <- match.arg(to)
  src_type <- attr(catalog, "type")

  if (to == "SBS96" && src_type == "SBS1536") {
    src_rows <- rownames(catalog)
    sbs96_rows <- paste0(substr(src_rows, 2L, 4L), substr(src_rows, 6L, 6L))
    out_rows <- mSigSpectra::catalog.row.order$SBS96
    out <- matrix(
      0, nrow = length(out_rows), ncol = ncol(catalog),
      dimnames = list(out_rows, colnames(catalog))
    )
    for (i in seq_along(sbs96_rows)) {
      out[sbs96_rows[i], ] <- out[sbs96_rows[i], ] + catalog[i, ]
    }
  } else if (to == "SBS96" && src_type == "SBS192") {
    # SBS192 rows are encoded with stranded trinucleotide context; collapse
    # by summing into their pyrimidine-canonical SBS96 rows.
    src_rows <- rownames(catalog)
    sbs96_rows <- pyr_tri(src_rows)
    out_rows <- mSigSpectra::catalog.row.order$SBS96
    out <- matrix(
      0, nrow = length(out_rows), ncol = ncol(catalog),
      dimnames = list(out_rows, colnames(catalog))
    )
    for (i in seq_along(sbs96_rows)) {
      out[sbs96_rows[i], ] <- out[sbs96_rows[i], ] + catalog[i, ]
    }
  } else if (to == "DBS78" && src_type == "DBS144") {
    src_rows <- rownames(catalog)
    dbs78_rows <- pyr_di(src_rows)
    out_rows <- mSigSpectra::catalog.row.order$DBS78
    out <- matrix(
      0, nrow = length(out_rows), ncol = ncol(catalog),
      dimnames = list(out_rows, colnames(catalog))
    )
    for (i in seq_along(dbs78_rows)) {
      out[dbs78_rows[i], ] <- out[dbs78_rows[i], ] + catalog[i, ]
    }
  } else {
    stop(
      "collapse_catalog: unsupported collapse ", src_type, " -> ", to
    )
  }

  as_catalog(
    out, type = to,
    ref_genome = attr(catalog, "ref_genome"),
    region     = attr(catalog, "region"),
    counts_or_density = attr(catalog, "counts_or_density")
  )
}

# SBS192 4-char row to SBS96 4-char row: if the central base is A or G,
# reverse-complement the trinucleotide and the ALT.
pyr_tri <- function(x) {
  stopifnot(all(nchar(x) == 4L))
  center <- substr(x, 2L, 2L)
  needs_rc <- center %in% c("A", "G")
  out <- x
  if (any(needs_rc)) {
    out[needs_rc] <- revc_sbs96(x[needs_rc])
  }
  out
}

# DBS144 4-char row to DBS78 4-char row: if the REF dinucleotide is not
# in canonical form, reverse-complement.
pyr_di <- function(x) {
  stopifnot(all(nchar(x) == 4L))
  # DBS78 canonical REFs from the shipped row order
  dbs78_refs <- unique(substr(mSigSpectra::catalog.row.order$DBS78, 1L, 2L))
  ref_di <- substr(x, 1L, 2L)
  needs_rc <- !(ref_di %in% dbs78_refs)
  out <- x
  if (any(needs_rc)) {
    out[needs_rc] <- revc_dbs144(x[needs_rc])
  }
  out
}
