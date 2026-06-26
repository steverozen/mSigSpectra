#!/usr/bin/env Rscript
# Build SBS96, DBS78, and ID83 catalogs from the duocarmycin Strelka VCFs
# in this folder, using mSigSpectra. Reference genome is GRCh37 (hs37d5).

suppressPackageStartupMessages({
  pkgload::load_all("/home/steve/github/mSigSpectra")
})

this_dir <- "/home/steve/github/mSigSpectra/development/duocarmycin-cell-line-exposures"

sbs_dbs_files <- list.files(
  this_dir, pattern = "_intersect\\.vcf\\.gz$", full.names = TRUE
)
sbs_dbs_files <- sbs_dbs_files[!grepl("INDELintersect", basename(sbs_dbs_files))]

id_files <- list.files(
  this_dir, pattern = "_INDELintersect\\.vcf\\.gz$", full.names = TRUE
)

sample_name_from <- function(file, suffix) {
  sub(paste0(suffix, "$"), "", basename(file))
}

process_sbs_dbs <- function(file) {
  sample <- sample_name_from(file, "_intersect.vcf.gz")
  message("=== SBS/DBS: ", sample)

  vcf <- read_vcf(file, filter = TRUE, name_of_vcf = sample)
  parts <- split_vcf(vcf, name_of_vcf = sample)

  if (nrow(parts$SBS) > 0L) {
    sbs_ann <- annotate_sbs_or_dbs_vcf(parts$SBS, ref_genome = "GRCh37")$annotated.vcf
    cat96 <- vcf_to_sbs_catalog(
      sbs_ann, type = "SBS96", ref_genome = "GRCh37",
      region = "genome", sample_name = sample
    )
    write_catalog(cat96, file.path(this_dir, paste0(sample, ".SBS96.csv")))
    message("  wrote SBS96 (", sum(cat96), " mutations)")
  } else {
    message("  no SBS variants")
  }

  if (nrow(parts$DBS) > 0L) {
    dbs_ann <- annotate_sbs_or_dbs_vcf(parts$DBS, ref_genome = "GRCh37")$annotated.vcf
    cat_dbs <- vcf_to_dbs_catalog(
      dbs_ann, type = "DBS78", ref_genome = "GRCh37",
      region = "genome", sample_name = sample
    )
    write_catalog(cat_dbs, file.path(this_dir, paste0(sample, ".DBS78.csv")))
    message("  wrote DBS78 (", sum(cat_dbs), " mutations)")
  } else {
    message("  no DBS variants")
  }
}

process_id <- function(file) {
  sample <- sample_name_from(file, "_INDELintersect.vcf.gz")
  message("=== ID: ", sample)

  vcf <- read_vcf(file, filter = TRUE, name_of_vcf = sample)
  parts <- split_vcf(vcf, name_of_vcf = sample)

  if (nrow(parts$ID) == 0L) {
    message("  no ID variants")
    return(invisible(NULL))
  }
  id_ann <- annotate_id_vcf(parts$ID, ref_genome = "GRCh37")$annotated.vcf
  for (type in c("ID83", "ID89", "ID476")) {
    cat_id <- vcf_to_id_catalog(
      id_ann, type = type, ref_genome = "GRCh37",
      region = "genome", sample_name = sample
    )
    out <- file.path(this_dir, paste0(sample, ".", type, ".csv"))
    if (type == "ID83") {
      write_catalog(cat_id, out)
    } else {
      # write_catalog has no row-header layout for ID89 / ID476 yet; emit a
      # minimal MutationType,<sample> CSV preserving canonical row order.
      df <- data.frame(MutationType = rownames(cat_id), cat_id, check.names = FALSE)
      utils::write.csv(df, out, row.names = FALSE)
    }
    message("  wrote ", type, " (", sum(cat_id), " mutations)")
  }
}

invisible(lapply(sbs_dbs_files, process_sbs_dbs))
invisible(lapply(id_files, process_id))

message("Done.")
