test_that("add_dna_region labels transcript variants G and others I", {
  vcf <- data.frame(
    CHROM = "1",
    POS = 1:4,
    trans.strand = c("+", "-", NA, "+")
  )
  out <- add_dna_region(vcf)
  expect_equal(out$dna.region, c("G", "G", "I", "G"))
  expect_equal(colnames(out)[ncol(out)], "dna.region")
})

test_that("add_dna_region leaves a VCF without trans.strand unchanged", {
  vcf <- data.frame(CHROM = "1", POS = 1:2)
  expect_message(out <- add_dna_region(vcf), "No trans.strand column")
  expect_false("dna.region" %in% colnames(out))
})

test_that("add_dna_region renames an existing dna.region column", {
  vcf <- data.frame(CHROM = "1", POS = 1:2, trans.strand = c("+", NA),
                    dna.region = c("x", "y"))
  expect_warning(out <- add_dna_region(vcf), "dna.region_old")
  expect_equal(out$dna.region_old, c("x", "y"))
  expect_equal(out$dna.region, c("G", "I"))
})

annotate_strelka_id <- function() {
  vcf <- read_vcf(extdata("Strelka-ID-GRCh37", "Strelka.ID.GRCh37.s1.vcf"),
                  filter = "PASS")
  sp <- suppressWarnings(split_vcf(vcf))
  suppressWarnings(annotate_id_vcf(sp$ID, ref_genome = "GRCh37"))
}

test_that("annotate_id_vcf adds dna.region only when the env var is TRUE", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))

  withr::local_envvar(MSIGSPECTRA_ADD_DNA_REGION = "TRUE")
  ann <- annotate_strelka_id()$annotated.vcf
  expect_true("dna.region" %in% colnames(ann))
  expect_setequal(unique(ann$dna.region), c("G", "I"))
  expect_equal(
    ann$dna.region == "G",
    ann$trans.strand %in% c("+", "-")
  )

  withr::local_envvar(MSIGSPECTRA_ADD_DNA_REGION = "FALSE")
  ann_off <- annotate_strelka_id()$annotated.vcf
  expect_false("dna.region" %in% colnames(ann_off))
  expect_equal(ann_off, ann[, setdiff(colnames(ann), "dna.region"),
                            with = FALSE])
})

test_that("annotate_id_vcf leaves dna.region off when the env var is unset", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))

  withr::local_envvar(MSIGSPECTRA_ADD_DNA_REGION = NA)
  ann <- annotate_strelka_id()$annotated.vcf
  expect_false("dna.region" %in% colnames(ann))
})
