#' Canonical row orders for catalogs
#'
#' Row-name orders for every catalog type supported by mSigSpectra. Exposed
#' for users who need to build catalogs from formats not otherwise supported.
#'
#' @format A named list of character vectors keyed by catalog type (e.g.
#'   `"SBS96"`, `"DBS78"`, `"ID83"`).
#'
#' @examples
#' catalog.row.order$SBS96[1:5]
#'
#' @name catalog.row.order
"catalog.row.order"

#' Transcript ranges for transcriptional strand annotation
#'
#' Precomputed transcript ranges (one row per gene) used by
#' [add_transcript_strand()] to determine the coding strand for each variant.
#' Sources: GENCODE v30 (human) and vM21 (mouse). Only genes with CCDS IDs are
#' retained.
#'
#' @format A [data.table::data.table] with columns `chrom`, `start`, `end`,
#'   `strand`, `Ensembl.gene.ID`, `gene.symbol`. One-based coordinates.
#'
#' @source <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gff3.gz>
#' @source <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gff3.gz>
#' @source <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz>
#'
#' @name trans.ranges
NULL

#' @rdname trans.ranges
"trans.ranges.GRCh37"

#' @rdname trans.ranges
"trans.ranges.GRCh38"

#' @rdname trans.ranges
"trans.ranges.GRCm38"

#' K-mer abundances for density calculations
#'
#' Nested list of k-mer counts keyed by
#' `BSgenome.Hsapiens.1000genomes.hs37d5` / `BSgenome.Hsapiens.UCSC.hg38` /
#' `BSgenome.Mmusculus.UCSC.mm10`, then by `exome` / `transcript` / `genome`,
#' then by catalog size as a string (`"78"`, `"96"`, `"136"`, `"144"`, `"192"`,
#' `"1536"`). Each leaf value is a named integer vector: k-mer → count.
#'
#' Used to convert counts to density (mutations per megabase of context).
#' See the planned `transform_catalog()` function (not yet ported).
#'
#' @examples
#' all.abundance$BSgenome.Hsapiens.UCSC.hg38$transcript$`144`[1:5]
#'
#' @name all.abundance
"all.abundance"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "all.abundance", "trans.ranges.GRCh37", "trans.ranges.GRCh38",
    "trans.ranges.GRCm38", "catalog.row.order",
    "CHROM", "POS", "POS2", "POS.plus.one", "POS.y",
    "REF", "REF.x", "REF.y", "ALT", "ALT.x", "ALT.y", "ref2alt",
    "FILTER", "ID", "SampleID",
    "chrom", "chrom.y", "strand", "bothstrand", "trans.strand",
    "trans.gene.symbol", "trans.Ensembl.gene.ID", "Ensembl.gene.ID",
    "gene.symbol", "trans.start.pos", "trans.end.pos",
    "exome.start", "exome.end", "dna.region",
    "seq.21bases", "pyr.mut", "mutation", "occurrences",
    "count", "n", "nrn", "N", "rn", "value", "variable",
    "type", "cols", "..col.names.order", "..column.to.use",
    "read.depth", "read.depth.x", "read.depth.y", "VAF", "VAF.x", "VAF.y",
    "delete.flag", "remark.for.DBS", "readthrough",
    "minus1bs", "minus2bs", "plus1bs", "plus2bs",
    "pos_id", "R", "HIGH", "LOW", "Exp_Level", "exp.value", "exp.level",
    "COSMIC_83", "Koh_89", "Koh_476",
    ".", ".orig_row"
  ))
}
