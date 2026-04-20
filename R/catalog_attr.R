# Mapping from number of rows to catalog type identifier.
.catalog_row_to_type <- c(
  "96"   = "SBS96",
  "192"  = "SBS192",
  "1536" = "SBS1536",
  "78"   = "DBS78",
  "136"  = "DBS136",
  "144"  = "DBS144",
  "83"   = "ID83",
  "166"  = "ID166",
  "1697" = "COMPOSITE"
)

# Legal catalog-type / counts-or-density values.
.legal_counts_or_density <- c(
  "counts", "density", "counts.signature", "density.signature"
)

# The `catalog.row.order` object uses key "ID" rather than "ID83"; map.
.catalog_row_order_key <- function(type) {
  if (type == "ID83") "ID" else type
}

#' Infer catalog type from the number of rows of a matrix
#'
#' @param n_rows Number of rows.
#'
#' @return A type identifier (e.g. `"SBS96"`, `"DBS78"`, `"ID83"`). Errors
#'   if `n_rows` does not correspond to a supported catalog type.
#'
#' @keywords internal
infer_catalog_type <- function(n_rows) {
  key <- as.character(n_rows)
  if (!key %in% names(.catalog_row_to_type)) {
    stop(
      "\nUnsupported number of rows for a catalog: ", n_rows,
      "\nExpected one of: ",
      paste(names(.catalog_row_to_type), collapse = ", ")
    )
  }
  .catalog_row_to_type[[key]]
}

#' Validate a counts_or_density argument
#'
#' @keywords internal
stop_if_counts_or_density_illegal <- function(counts_or_density) {
  if (is.null(counts_or_density)) return(invisible(NULL))
  if (!counts_or_density %in% .legal_counts_or_density) {
    stop(
      "Unrecognized counts_or_density identifier: ", counts_or_density,
      "\nNeed one of ",
      paste(.legal_counts_or_density, collapse = ", ")
    )
  }
  invisible(NULL)
}

#' Reorder catalog rows to the canonical order for its type
#'
#' @param x A matrix with rownames.
#' @param type Catalog type (e.g. `"SBS96"`).
#'
#' @return `x` with rows reordered to match `catalog.row.order[[type]]`.
#'   Errors if any canonical rowname is missing.
#'
#' @keywords internal
check_and_reorder_rownames <- function(x, type) {
  key <- .catalog_row_order_key(type)
  correct <- mSigSpectra::catalog.row.order[[key]]
  missing_rows <- setdiff(correct, rownames(x))
  if (length(missing_rows) > 0L) {
    stop(
      "\nThe input matrix does not have correct rownames for a ", type,
      " catalog.\nMissing rownames:\n",
      paste(missing_rows, collapse = " ")
    )
  }
  x[correct, , drop = FALSE]
}

#' Infer k-mer abundance from attributes
#'
#' Looks up the appropriate entry of `all.abundance` keyed on
#' `ref_genome` and `region`. Returns `NULL` if no entry is found,
#' if the catalog is a signature / density catalog with a flat abundance,
#' or if the catalog is COMPOSITE (which has no meaningful abundance).
#'
#' @keywords internal
infer_abundance <- function(x, ref_genome, region, counts_or_density) {
  stop_if_region_illegal(region)
  stop_if_counts_or_density_illegal(counts_or_density)

  if (nrow(x) == 1697L) return(NULL)

  if (is_density(counts_or_density)) {
    # Flat abundance (one per context size)
    flat <- get("flat.abundance", envir = asNamespace("mSigSpectra"))
    ab <- flat[[as.character(nrow(x))]]
    return(ab)
  }

  if (is.null(ref_genome)) return(NULL)

  ref_genome_name <- if (inherits(ref_genome, "BSgenome")) {
    ref_genome@pkgname
  } else {
    infer_ref_genome_name(ref_genome)
  }

  ab <- mSigSpectra::all.abundance[[ref_genome_name]]
  if (is.null(ab)) return(NULL)

  ab_region <- ab[[region]]
  if (is.null(ab_region)) return(NULL)

  ab_region[[as.character(nrow(x))]]
}

is_density <- function(counts_or_density) {
  !is.null(counts_or_density) &&
    counts_or_density %in% c("density", "density.signature")
}

#' Map a character ref_genome argument to its canonical BSgenome package name
#'
#' @keywords internal
infer_ref_genome_name <- function(ref_genome) {
  if (is.null(ref_genome)) stop("Need a non-NULL ref_genome")
  if (inherits(ref_genome, "BSgenome")) return(ref_genome@pkgname)

  stopifnot(is.character(ref_genome), length(ref_genome) == 1L)

  if (ref_genome %in% c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
    "BSgenome.Hsapiens.UCSC.hg38"
  } else if (ref_genome %in%
             c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
    "BSgenome.Hsapiens.1000genomes.hs37d5"
  } else if (ref_genome %in%
             c("GRCm38", "mm10", "BSgenome.Mmusculus.UCSC.mm10")) {
    "BSgenome.Mmusculus.UCSC.mm10"
  } else {
    stop(
      "Unrecognized ref_genome: ", ref_genome,
      "\nNeed one of GRCh38/hg38/GRCh37/hg19/GRCm38/mm10"
    )
  }
}

#' Turn a numeric matrix into a mutational-spectrum catalog
#'
#' Attaches the standard catalog attributes (`type`, `counts_or_density`,
#' `ref_genome`, `region`, `abundance`) to `x` and returns it. **No S3
#' class is set** — mSigSpectra catalogs are plain matrices with
#' attributes; functions key off `attr(x, "type")` via explicit checks.
#'
#' @param x A numeric matrix, data.frame coercible to numeric, or a named
#'   numeric vector (converted to a one-column matrix with the vector names
#'   as rownames).
#' @param type Optional catalog type identifier (`"SBS96"`, `"DBS78"`,
#'   `"ID83"`, etc.). Inferred from `nrow(x)` if `NULL`.
#' @param ref_genome Optional BSgenome object or alias; recorded as an
#'   attribute and used for abundance lookup.
#' @param region One of `"genome"`, `"exome"`, `"transcript"`,
#'   `"unknown"`. For stranded catalogs (`SBS192`, `DBS144`) `"genome"` is
#'   silently promoted to `"transcript"`.
#' @param abundance Optional named numeric vector of k-mer counts. If
#'   `NULL`, inferred from shipped `all.abundance` keyed on `ref_genome`
#'   and `region`.
#' @param counts_or_density One of `"counts"`, `"density"`,
#'   `"counts.signature"`, `"density.signature"`.
#' @param infer_rownames If `TRUE` and `x` has no rownames, the canonical
#'   rownames for the catalog type are attached (**assuming the row order
#'   is already correct**). If `FALSE`, `x` must already have the canonical
#'   rownames.
#'
#' @return `x` as a numeric matrix with attributes set.
#'
#' @examples
#' m <- matrix(
#'   1, nrow = 96, ncol = 1,
#'   dimnames = list(catalog.row.order$SBS96, "sample1")
#' )
#' cat96 <- as_catalog(m)
#' attr(cat96, "type")   # "SBS96"
#' attr(cat96, "region") # "unknown"
#'
#' @export
as_catalog <- function(x,
                       type = NULL,
                       ref_genome = NULL,
                       region = "unknown",
                       abundance = NULL,
                       counts_or_density = "counts",
                       infer_rownames = FALSE) {

  if (!is.matrix(x)) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    } else if (is.vector(x)) {
      nms <- names(x)
      x <- matrix(x, ncol = 1L)
      rownames(x) <- nms
      colnames(x) <- "Unknown"
    } else {
      stop("x must be a numeric matrix, data.frame, or named numeric vector")
    }
  }
  storage.mode(x) <- "numeric"

  if (is.null(type)) {
    type <- infer_catalog_type(nrow(x))
  } else if (!type %in% .catalog_row_to_type) {
    stop(
      "Unrecognized catalog type: ", type,
      "\nNeed one of ",
      paste(unique(.catalog_row_to_type), collapse = ", ")
    )
  }

  if (is.null(rownames(x))) {
    if (!infer_rownames) {
      stop(
        "x has no rownames. Either set the canonical rownames first, ",
        "or pass infer_rownames = TRUE and ensure row order is correct."
      )
    }
    rownames(x) <- mSigSpectra::catalog.row.order[[.catalog_row_order_key(type)]]
  } else {
    x <- check_and_reorder_rownames(x, type)
  }

  stop_if_region_illegal(region)
  stop_if_counts_or_density_illegal(counts_or_density)

  if (type %in% c("SBS192", "DBS144")) {
    if (region == "genome") region <- "transcript"
    stop_if_transcribed_region_illegal(region)
  }

  if (!is.null(ref_genome)) {
    ref_genome <- normalize_genome_arg(ref_genome)
  }

  if (is.null(abundance)) {
    abundance <- infer_abundance(x, ref_genome, region, counts_or_density)
  }

  attr(x, "type")              <- type
  attr(x, "counts_or_density") <- counts_or_density
  attr(x, "ref_genome")        <- ref_genome
  attr(x, "region")            <- region
  attr(x, "abundance")         <- abundance

  x
}

#' Check whether an object looks like an mSigSpectra catalog
#'
#' Returns `TRUE` if `x` is a numeric matrix with the five catalog
#' attributes (`type`, `counts_or_density`, `ref_genome`, `region`,
#' `abundance`) and canonical rownames for its type.
#'
#' @param x Any R object.
#'
#' @export
is_catalog <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) return(FALSE)
  required_attrs <- c("type", "counts_or_density", "region")
  if (!all(required_attrs %in% names(attributes(x)))) return(FALSE)
  type <- attr(x, "type")
  if (!type %in% .catalog_row_to_type) return(FALSE)
  if (nrow(x) != as.integer(names(.catalog_row_to_type)[.catalog_row_to_type == type])) {
    return(FALSE)
  }
  key <- .catalog_row_order_key(type)
  correct <- mSigSpectra::catalog.row.order[[key]]
  identical(rownames(x), correct)
}

#' Report the attributes of an mSigSpectra catalog
#'
#' @param x A catalog.
#'
#' @return A named list with elements `type`, `counts_or_density`,
#'   `ref_genome`, `region`, `abundance`.
#'
#' @export
catalog_attrs <- function(x) {
  list(
    type              = attr(x, "type"),
    counts_or_density = attr(x, "counts_or_density"),
    ref_genome        = attr(x, "ref_genome"),
    region            = attr(x, "region"),
    abundance         = attr(x, "abundance")
  )
}

#' Combine catalogs across samples (column-bind)
#'
#' Checks that all input catalogs share the same `type`,
#' `counts_or_density`, `ref_genome`, and `region`, then `cbind`s their
#' matrices and re-applies the shared attributes.
#'
#' @param catalogs A list of catalogs.
#'
#' @return A single catalog with `ncol` equal to the sum of the input
#'   `ncol`s.
#'
#' @export
cbind_catalogs <- function(catalogs) {
  stopifnot(is.list(catalogs), length(catalogs) >= 1L)

  ref <- catalog_attrs(catalogs[[1L]])
  for (i in seq_along(catalogs)[-1L]) {
    this <- catalog_attrs(catalogs[[i]])
    if (!identical(this$type, ref$type) ||
        !identical(this$counts_or_density, ref$counts_or_density) ||
        !identical(this$region, ref$region)) {
      stop(
        "Cannot cbind catalogs with mismatched ",
        "type / counts_or_density / region attributes"
      )
    }
  }

  mat <- do.call(cbind, lapply(catalogs, function(x) {
    storage.mode(x) <- "numeric"
    m <- unclass(x)
    attributes(m) <- list(dim = dim(m), dimnames = dimnames(m))
    m
  }))

  attr(mat, "type")              <- ref$type
  attr(mat, "counts_or_density") <- ref$counts_or_density
  attr(mat, "ref_genome")        <- ref$ref_genome
  attr(mat, "region")            <- ref$region
  attr(mat, "abundance")         <- ref$abundance
  mat
}

#' Subset a catalog while preserving attributes
#'
#' Base `[` on a matrix drops attributes other than `dim` / `dimnames`.
#' Use `subset_catalog()` when you want to keep the catalog's type /
#' ref_genome / region / abundance metadata through a subset operation.
#'
#' @param x A catalog.
#' @param rows,cols Numeric, logical, or character indices (see
#'   [base::Extract]). If `NULL`, all rows / columns are kept.
#'
#' @export
subset_catalog <- function(x, rows = NULL, cols = NULL) {
  attrs <- catalog_attrs(x)
  if (is.null(rows)) rows <- seq_len(nrow(x))
  if (is.null(cols)) cols <- seq_len(ncol(x))
  m <- unclass(x)[rows, cols, drop = FALSE]
  # If rows were subset, abundance may no longer be meaningful
  if (length(rows) != nrow(x)) attrs$abundance <- NULL
  attr(m, "type")              <- attrs$type
  attr(m, "counts_or_density") <- attrs$counts_or_density
  attr(m, "ref_genome")        <- attrs$ref_genome
  attr(m, "region")            <- attrs$region
  attr(m, "abundance")         <- attrs$abundance
  m
}
