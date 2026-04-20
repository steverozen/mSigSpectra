#' Normalize a reference-genome argument to a BSgenome object
#'
#' Accepts either a `BSgenome` object directly, or one of the string
#' identifiers `"GRCh37"` / `"hg19"` / `"BSgenome.Hsapiens.1000genomes.hs37d5"`,
#' `"GRCh38"` / `"hg38"` / `"BSgenome.Hsapiens.UCSC.hg38"`, or
#' `"GRCm38"` / `"mm10"` / `"BSgenome.Mmusculus.UCSC.mm10"`. The relevant
#' BSgenome package must be installed; an informative error is raised if not.
#'
#' @param ref_genome A BSgenome object or a recognized character identifier.
#'
#' @return A BSgenome object.
#'
#' @keywords internal
normalize_genome_arg <- function(ref_genome) {
  if (is.null(ref_genome)) stop("Need a non-NULL ref_genome")
  if (inherits(ref_genome, "BSgenome")) return(ref_genome)

  stopifnot(is.character(ref_genome), length(ref_genome) == 1L)

  if (ref_genome %in%
      c("GRCh38", "hg38", "BSgenome.Hsapiens.UCSC.hg38")) {
    pkg <- "BSgenome.Hsapiens.UCSC.hg38"
  } else if (ref_genome %in%
             c("GRCh37", "hg19", "BSgenome.Hsapiens.1000genomes.hs37d5")) {
    pkg <- "BSgenome.Hsapiens.1000genomes.hs37d5"
  } else if (ref_genome %in%
             c("GRCm38", "mm10", "BSgenome.Mmusculus.UCSC.mm10")) {
    pkg <- "BSgenome.Mmusculus.UCSC.mm10"
  } else {
    stop(
      "Unrecognized ref_genome: ", ref_genome,
      "\nNeed a BSgenome reference genome or one of the strings ",
      "GRCh38, hg38, GRCh37, hg19, GRCm38, mm10"
    )
  }

  if ("" == system.file(package = pkg)) {
    stop(
      "\nPlease install ", pkg, ":\n",
      "BiocManager::install(\"", pkg, "\")"
    )
  }
  stopifnot(requireNamespace(pkg, quietly = TRUE))
  getExportedValue(pkg, pkg)
}

#' Validate a `region` argument
#'
#' @param region Character string to check.
#'
#' @return `NULL` invisibly; raises an error if `region` is not one of
#'   `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.
#'
#' @keywords internal
stop_if_region_illegal <- function(region) {
  if (!region %in% c("genome", "exome", "transcript", "unknown")) {
    stop(
      "Unrecognized region identifier: ", region,
      "\nNeed one of genome, exome, transcript, unknown"
    )
  }
  invisible(NULL)
}

#' Validate a `region` argument for catalog types that require transcript strand
#'
#' SBS192, DBS144 and similar stranded catalogs cannot be built from
#' `region = "genome"`, since variants outside transcripts have no strand.
#'
#' @param region Character string to check.
#'
#' @keywords internal
stop_if_transcribed_region_illegal <- function(region) {
  if (!region %in% c("transcript", "exome", "unknown")) {
    stop(
      "Require region to be one of transcript, exome, or unknown ",
      "for stranded catalogs (SBS192, DBS144). Got ", region
    )
  }
  invisible(NULL)
}
