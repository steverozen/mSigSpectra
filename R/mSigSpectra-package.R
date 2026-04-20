#' mSigSpectra: Build mutational-spectrum catalogs from VCF files
#'
#' Reads variant call files (VCFs) caller-agnostically, annotates variants with
#' flanking sequence context and transcriptional strand, and builds mutational-
#' spectrum catalogs (SBS, DBS, indel) at multiple resolutions and in both
#' counts and density representations.
#'
#' Catalogs are plain numeric matrices (rows = mutation categories,
#' columns = samples) carrying attributes `type`, `ref_genome`, `region`,
#' `counts_or_density`, and optionally `abundance`. mSigSpectra deliberately
#' does not use S3 classes or method dispatch on catalogs; behavior is selected
#' by explicit `attr(x, "type")` checks inside regular functions.
#'
#' Plotting is intentionally out of scope and is handled by the separate
#' mSigPlot package.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib mSigSpectra, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
