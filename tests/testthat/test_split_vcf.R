test_that("split_vcf classifies rows by REF/ALT length alone", {
  vcf <- data.frame(
    CHROM = rep("1", 6),
    POS   = 1:6,
    REF   = c("A", "AT", "ATG", "A", "ATC", "A"),
    ALT   = c("G", "GC", "CCC", "ATG", "AT", "GC")
  )
  out <- suppressWarnings(split_vcf(vcf))
  expect_equal(nrow(out$SBS), 1L)           # A>G
  expect_equal(nrow(out$DBS), 1L)           # AT>GC
  expect_equal(nrow(out$ID),  3L)           # A>ATG, ATC>AT, A>GC (all unequal)
  expect_equal(nrow(out$discarded), 1L)     # ATG>CCC (3/3 substitution)
  expect_true("discarded.reason" %in% colnames(out$discarded))
})

test_that("split_vcf returns empty tables + NULL discarded for an empty VCF", {
  vcf <- data.frame(CHROM = character(), POS = integer(),
                    REF = character(), ALT = character())
  out <- split_vcf(vcf)
  expect_equal(nrow(out$SBS), 0L)
  expect_equal(nrow(out$DBS), 0L)
  expect_equal(nrow(out$ID),  0L)
  expect_null(out$discarded)
})

test_that("split_vcf returns NULL discarded when no rows are discarded", {
  vcf <- data.frame(CHROM = "1", POS = 1L, REF = "A", ALT = "G")
  out <- split_vcf(vcf)
  expect_equal(nrow(out$SBS), 1L)
  expect_null(out$discarded)
})

test_that("split_vcf works on a real Strelka ID VCF (all rows classified as ID)", {
  vcf <- read_vcf("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")
  out <- suppressWarnings(split_vcf(vcf))
  expect_equal(nrow(out$SBS), 0L)
  expect_equal(nrow(out$DBS), 0L)
  expect_gt(nrow(out$ID),  0L)
  # All rows have unequal REF/ALT lengths
  expect_true(all(nchar(out$ID$REF) != nchar(out$ID$ALT)))
})

test_that("split_vcf works on a real Strelka SBS VCF (all rows classified as SBS)", {
  vcf <- read_vcf("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf")
  out <- suppressWarnings(split_vcf(vcf))
  expect_gt(nrow(out$SBS), 0L)
  expect_equal(nrow(out$ID), 0L)
})
