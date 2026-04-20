# Expected first-column header for each ICAMS-native catalog type.
# Each entry is the number of leading "header" columns the CSV file has
# (before the sample columns) and their expected column names.
.icams_native_headers <- list(
  SBS96    = list(n = 2L, cols = c("Mutation type", "Trinucleotide")),
  SBS192   = list(n = 3L, cols = c("Strand",        "Mutation type", "Trinucleotide")),
  SBS1536  = list(n = 2L, cols = c("Mutation type", "Pentanucleotide")),
  DBS78    = list(n = 2L, cols = c("Ref",           "Var")),
  DBS136   = list(n = 1L, cols = c("Quad")),
  DBS144   = list(n = 2L, cols = c("Ref",           "Var")),
  ID83     = list(n = 4L, cols = c("Type", "Subtype", "Indel_size", "Repeat_MH_size")),
  ID166    = list(n = 5L, cols = c("Region", "Type", "Subtype", "Indel_size", "Repeat_MH_size"))
)

#' Read a mutational-spectrum catalog from a file
#'
#' Auto-detects the file format from its first row. Recognized formats:
#' * **ICAMS native** CSV — first N columns are header columns (mutation type,
#'   context, etc), remaining columns are per-sample counts / densities.
#' * **SigProfiler** TSV / CSV — single header column with bracketed mutation
#'   labels like `A[C>A]A`, `AA[C>A]AA`, or PCAWG indel codes like
#'   `1:Del:C:0`.
#' * **COSMIC** CSV — SBS96 uses the ICAMS-external row-header format;
#'   SBS192 uses a stranded variant of the same; DBS78 / ID83 also share
#'   their SigProfiler-style layout.
#'
#' The matrix is reordered to the canonical rownames for its type and
#' wrapped in a catalog via [as_catalog()].
#'
#' @param file Path to the catalog file.
#' @param ref_genome Optional BSgenome object or alias; stored as attribute.
#' @param region One of `"genome"`, `"exome"`, `"transcript"`,
#'   `"unknown"`.
#' @param counts_or_density One of `"counts"`, `"density"`,
#'   `"counts.signature"`, `"density.signature"`.
#' @param format `"auto"` (default), `"ICAMS"`, `"SigProfiler"`, or
#'   `"COSMIC"`. In `"auto"` mode the format is inferred from the file's
#'   first data row.
#'
#' @return A catalog matrix with attributes (see [as_catalog()]).
#'
#' @export
read_catalog <- function(file,
                         ref_genome = NULL,
                         region = "unknown",
                         counts_or_density = "counts",
                         format = c("auto", "ICAMS", "SigProfiler", "COSMIC")) {
  format <- match.arg(format)
  stop_if_region_illegal(region)
  stop_if_counts_or_density_illegal(counts_or_density)

  dt <- data.table::fread(file)
  # Drop all-NA columns (some writers emit trailing empty columns).
  keep_cols <- which(vapply(dt, function(x) !all(is.na(x)), logical(1)))
  dt <- dt[, keep_cols, with = FALSE]
  # Drop rows where all values are NA.
  dt <- stats::na.omit(dt)

  if (format == "auto") {
    format <- detect_catalog_format(dt)
  }

  m <- switch(
    format,
    ICAMS       = read_icams_native(dt),
    SigProfiler = read_sigprofiler(dt),
    COSMIC      = read_cosmic(dt)
  )

  as_catalog(
    m,
    ref_genome = ref_genome,
    region = region,
    counts_or_density = counts_or_density
  )
}

detect_catalog_format <- function(dt) {
  first_col <- as.character(unlist(dt[, 1L]))
  # SigProfiler-style: "A[C>A]A", "AA[C>A]AA", or PCAWG indel codes.
  if (any(grepl("\\[", first_col))) return("SigProfiler")
  if (any(c("1:Del:C:0", "1:Del:T:0") %in% first_col)) return("SigProfiler")
  # Otherwise assume ICAMS-native / COSMIC CSV (both use the same layout).
  "ICAMS"
}

read_icams_native <- function(dt) {
  n_rows <- nrow(dt)
  type <- infer_catalog_type(n_rows)
  header <- .icams_native_headers[[type]]
  if (is.null(header)) {
    stop("read_catalog: no native header layout for type ", type)
  }

  n_header_cols <- header$n
  data_cols <- seq.int(n_header_cols + 1L, ncol(dt))
  if (length(data_cols) < 1L) {
    stop("read_catalog: catalog has no sample columns")
  }

  m <- as.matrix(dt[, data_cols, with = FALSE])
  rns <- build_icams_native_rownames(dt[, seq_len(n_header_cols), with = FALSE],
                                     type)
  rownames(m) <- rns

  # Coerce to numeric matrix
  storage.mode(m) <- "numeric"

  # Canonicalize row order
  key <- if (type == "ID83") "ID" else type
  correct <- mSigSpectra::catalog.row.order[[key]]
  if (!setequal(rownames(m), correct)) {
    stop(
      "read_catalog: catalog rownames do not match canonical order for ",
      type
    )
  }
  m[correct, , drop = FALSE]
}

build_icams_native_rownames <- function(header_dt, type) {
  # SBS96: row name = paste0(trinucleotide, ALT), where "Mutation type" is
  # e.g. "C>A" and "Trinucleotide" is e.g. "ACA" -> row "ACAA".
  if (type == "SBS96") {
    mut <- as.character(unlist(header_dt[[1L]]))
    tri <- as.character(unlist(header_dt[[2L]]))
    alt <- substr(mut, 3L, 3L)
    return(paste0(tri, alt))
  }
  if (type == "SBS1536") {
    mut <- as.character(unlist(header_dt[[1L]]))
    pen <- as.character(unlist(header_dt[[2L]]))
    alt <- substr(mut, 3L, 3L)
    return(paste0(pen, alt))
  }
  if (type == "SBS192") {
    strand <- as.character(unlist(header_dt[[1L]]))
    mut    <- as.character(unlist(header_dt[[2L]]))
    tri    <- as.character(unlist(header_dt[[3L]]))
    alt    <- substr(mut, 3L, 3L)
    # For transcribed strand ("T"), reverse-complement trinucleotide + alt
    rows <- paste0(tri, alt)
    rc_idx <- which(strand == "T")
    if (length(rc_idx) > 0L) {
      rows[rc_idx] <- revc_sbs96(rows[rc_idx])
    }
    return(rows)
  }
  if (type == "DBS78" || type == "DBS144") {
    ref <- as.character(unlist(header_dt[[1L]]))
    var <- as.character(unlist(header_dt[[2L]]))
    return(paste0(ref, var))
  }
  if (type == "DBS136") {
    return(as.character(unlist(header_dt[[1L]])))
  }
  if (type == "ID83") {
    # 4 header columns: Type, Subtype, Indel_size, Repeat_MH_size
    # Row name format is e.g. "DEL:T:1:2" -> concatenate with ":"
    a <- as.character(unlist(header_dt[[1L]]))
    b <- as.character(unlist(header_dt[[2L]]))
    c <- as.character(unlist(header_dt[[3L]]))
    d <- as.character(unlist(header_dt[[4L]]))
    return(paste(a, b, c, d, sep = ":"))
  }
  if (type == "ID166") {
    # Header columns: Region, Type, Subtype, Indel_size, Repeat_MH_size.
    # Row labels in catalog.row.order$ID166 start with the region letter
    # (e.g. "G" for genic) followed by ":" and the ID83-style key, so
    # the row label is paste(Region, Type, Subtype, Indel_size, Repeat_MH_size).
    r <- as.character(unlist(header_dt[[1L]]))
    a <- as.character(unlist(header_dt[[2L]]))
    b <- as.character(unlist(header_dt[[3L]]))
    c <- as.character(unlist(header_dt[[4L]]))
    d <- as.character(unlist(header_dt[[5L]]))
    return(paste(r, a, b, c, d, sep = ":"))
  }
  stop("build_icams_native_rownames: no layout for type ", type)
}

# SigProfiler reader (SBS96 bracketed form: "A[C>A]A"; SBS1536
# "AA[C>A]AA"; ID83 PCAWG form: "1:Del:C:0"). COSMIC reader is a thin
# pass-through to read_icams_native for the SBS96/SBS192/DBS78 cases
# (they share the row-header layout).

read_sigprofiler <- function(dt) {
  first_col <- as.character(unlist(dt[, 1L]))
  if (any(grepl("^[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]$", first_col))) {
    return(read_sigprofiler_sbs96(dt))
  }
  if (any(grepl("^[ACGT]{2}\\[[ACGT]>[ACGT]\\][ACGT]{2}$", first_col))) {
    return(read_sigprofiler_sbs1536(dt))
  }
  if (any(grepl(":(Del|Ins):", first_col))) {
    return(read_sigprofiler_id83(dt))
  }
  stop("read_sigprofiler: unrecognized first-column format")
}

read_sigprofiler_sbs96 <- function(dt) {
  labels <- as.character(unlist(dt[, 1L]))
  # "A[C>A]A" -> context "ACA" alt "A" -> "ACAA"
  trin <- paste0(substr(labels, 1L, 1L), substr(labels, 3L, 3L),
                 substr(labels, 7L, 7L))
  alt  <- substr(labels, 5L, 5L)
  rns  <- paste0(trin, alt)
  m <- as.matrix(dt[, -1L, with = FALSE])
  rownames(m) <- rns
  storage.mode(m) <- "numeric"
  if (!setequal(rownames(m), mSigSpectra::catalog.row.order$SBS96)) {
    stop("read_sigprofiler_sbs96: row labels do not cover SBS96")
  }
  m[mSigSpectra::catalog.row.order$SBS96, , drop = FALSE]
}

read_sigprofiler_sbs1536 <- function(dt) {
  labels <- as.character(unlist(dt[, 1L]))
  # "AA[C>A]AA" -> pent "AACAA" alt "A"
  pent <- paste0(substr(labels, 1L, 2L),
                 substr(labels, 4L, 4L),
                 substr(labels, 8L, 9L))
  alt  <- substr(labels, 6L, 6L)
  rns  <- paste0(pent, alt)
  m <- as.matrix(dt[, -1L, with = FALSE])
  rownames(m) <- rns
  storage.mode(m) <- "numeric"
  if (!setequal(rownames(m), mSigSpectra::catalog.row.order$SBS1536)) {
    stop("read_sigprofiler_sbs1536: row labels do not cover SBS1536")
  }
  m[mSigSpectra::catalog.row.order$SBS1536, , drop = FALSE]
}

read_sigprofiler_id83 <- function(dt) {
  # ICAMS ships a SigPro.to.ICAMS.ID mapping in sysdata.rda. Use it if the
  # incoming labels cover the 96-row "ID96" SigProfiler layout; otherwise
  # assume labels are already in ICAMS ID83 form.
  labels <- as.character(unlist(dt[, 1L]))

  if (all(labels %in% mSigSpectra::catalog.row.order$ID)) {
    m <- as.matrix(dt[, -1L, with = FALSE])
    rownames(m) <- labels
    storage.mode(m) <- "numeric"
    return(m[mSigSpectra::catalog.row.order$ID, , drop = FALSE])
  }

  # Drop labels that are not part of ID83 (complex, non_matching, Ins:M*)
  drop <- grepl("^(complex|non[._-]?matching)$", labels, ignore.case = TRUE) |
    grepl(":Ins:M", labels)
  if (any(drop)) {
    dt <- dt[!drop, ]
    labels <- labels[!drop]
  }

  mapping <- tryCatch(
    get("SigPro.to.ICAMS.ID", envir = asNamespace("mSigSpectra")),
    error = function(e) NULL
  )
  if (is.null(mapping)) {
    stop(
      "read_sigprofiler_id83: cannot translate SigProfiler ID labels ",
      "-- SigPro.to.ICAMS.ID mapping not found in sysdata"
    )
  }

  icams_labels <- mapping[labels]
  if (any(is.na(icams_labels))) {
    stop(
      "read_sigprofiler_id83: unrecognized SigProfiler ID labels: ",
      paste(labels[is.na(icams_labels)], collapse = ", ")
    )
  }

  m <- as.matrix(dt[, -1L, with = FALSE])
  rownames(m) <- icams_labels
  storage.mode(m) <- "numeric"
  if (!setequal(rownames(m), mSigSpectra::catalog.row.order$ID)) {
    stop("read_sigprofiler_id83: row labels do not cover ID83")
  }
  m[mSigSpectra::catalog.row.order$ID, , drop = FALSE]
}

read_cosmic <- function(dt) {
  # COSMIC SBS96 / SBS192 / DBS78 share the ICAMS-external layout.
  read_icams_native(dt)
}

#' Write a mutational-spectrum catalog to a file
#'
#' Writes `catalog` in ICAMS-native CSV format (the only format currently
#' supported for writing). The row-header columns that precede the sample
#' columns are taken from the shipped `catalog.row.headers` object for the
#' catalog's type.
#'
#' @param catalog A catalog (see [as_catalog()]).
#' @param file Output path.
#' @param format Output format. Currently only `"ICAMS"` is supported;
#'   `"SigProfiler"` and `"COSMIC"` error with an informative message.
#' @param sep Column separator. Defaults to `","` for ICAMS format.
#'
#' @return `file` invisibly.
#'
#' @export
write_catalog <- function(catalog, file,
                          format = c("ICAMS", "SigProfiler", "COSMIC"),
                          sep = ",") {
  format <- match.arg(format)
  if (format != "ICAMS") {
    stop(
      "write_catalog(format = '", format, "') is not yet implemented. ",
      "Currently only format = 'ICAMS' is supported."
    )
  }

  type <- attr(catalog, "type")
  if (is.null(type)) {
    stop(
      "catalog has no 'type' attribute; use as_catalog() to attach one ",
      "before writing."
    )
  }
  key <- if (type == "ID83") "ID" else type

  # Look up the shipped header rows for this type.
  headers <- get("catalog.row.headers", envir = asNamespace("mSigSpectra"))
  row_header <- headers[[key]]
  if (is.null(row_header)) {
    stop("write_catalog: no row header layout for type ", type)
  }

  # Reorder catalog to canonical row order.
  correct <- mSigSpectra::catalog.row.order[[key]]
  if (!setequal(rownames(catalog), correct)) {
    stop("write_catalog: catalog rownames do not match canonical ", type)
  }
  catalog_ordered <- catalog[correct, , drop = FALSE]

  # The shipped header rows are already in canonical row order; cbind them
  # alongside the numeric sample columns and fwrite.
  data_dt <- data.table::as.data.table(catalog_ordered)
  combined <- cbind(data.table::as.data.table(row_header), data_dt)
  data.table::fwrite(combined, file = file, sep = sep)
  invisible(file)
}
