#!/usr/bin/env Rscript
# Create inst/extdata/annotated-ID-GRCh37/HepG2_Duo_100pM_2mth_cl2.annotated.indel.vcf.gz,
# the pre-annotated indel VCF used by tests/testthat/test_build_one_file_rosetta.R,
# from the duocarmycin cell-line Strelka indel VCF. Needs the hs37d5 BSgenome.
devtools::load_all(quiet = TRUE)
f <- here::here(
  "development", "duocarmycin-cell-line-exposures",
  "HepG2_Duo_100pM_2mth_cl2_INDELintersect.vcf.gz"
)
vcf <- read_vcf(f, filter = TRUE, name_of_vcf = "HepG2_cl2")
sp <- split_vcf(vcf, name_of_vcf = "HepG2_cl2")
ann <- annotate_id_vcf(sp$ID, ref_genome = "GRCh37")$annotated.vcf
out <- here::here(
  "inst", "extdata", "annotated-ID-GRCh37",
  "HepG2_Duo_100pM_2mth_cl2.annotated.indel.vcf.gz"
)
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(out, "w")
utils::write.table(ann, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
message("Wrote ", out, " (", nrow(ann), " indels)")
