test_that("standard_chrom_name strips chr prefix and drops non-canonical chromosomes", {
  df <- data.frame(
    chrom = c("chr1", "2", "GL000001", "KI270711", "chrM", "chrX", "JH001"),
    pos = 1:7
  )
  out <- standard_chrom_name(df)
  expect_equal(out$chrom, c("1", "2", "X"))
})

test_that("standard_chrom_name_new returns a list with df + discarded.variants", {
  df <- data.frame(
    CHROM = c("1", "2", "GL000001", "KI270711", "random1", "Hs37d5", "chrM", "JH001"),
    POS = 1:8
  )
  expect_warning(out <- standard_chrom_name_new(df), "GL")
  expect_equal(nrow(out$df), 2)  # "1" and "2"
  expect_equal(nrow(out$discarded.variants), 6)
  expect_true("discarded.reason" %in% names(out$discarded.variants))
})

test_that("standard_chrom_name_new returns no discarded element when all rows are standard", {
  df <- data.frame(CHROM = c("1", "2", "X", "Y"), POS = 1:4)
  out <- standard_chrom_name_new(df)
  expect_equal(nrow(out$df), 4)
  expect_null(out$discarded.variants)
})

test_that("select_variants_by_chrom_name filters rows and reports discards", {
  df <- data.frame(CHROM = c("1", "2", "3", "1", "X"), POS = 1:5)
  out <- suppressWarnings(select_variants_by_chrom_name(df, c("1", "2")))
  expect_equal(nrow(out$df), 3)
  expect_equal(nrow(out$discarded.variants), 2)
  expect_equal(
    unique(out$discarded.variants$discarded.reason),
    "Chromosome names not selected by user"
  )
})

test_that("select_variants_by_chrom_name returns no discards when all rows pass", {
  df <- data.frame(CHROM = c("1", "2"), POS = 1:2)
  out <- select_variants_by_chrom_name(df, c("1", "2", "3"))
  expect_equal(nrow(out$df), 2)
  expect_null(out$discarded.variants)
})

test_that("check_and_fix_chrom_names adds chr prefix when ref has prefix and VCF does not", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  ref <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  vcf.df <- data.frame(CHROM = c("1", "2", "X"), POS = 1:3)
  out <- check_and_fix_chrom_names(vcf.df, ref)
  expect_equal(out, c("chr1", "chr2", "chrX"))
})

test_that("check_and_fix_chrom_names strips chr prefix when ref lacks it", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  ref <- BSgenome::getBSgenome("BSgenome.Hsapiens.1000genomes.hs37d5")
  vcf.df <- data.frame(CHROM = c("chr1", "chr2", "chrX"), POS = 1:3)
  out <- check_and_fix_chrom_names(vcf.df, ref)
  expect_equal(out, c("1", "2", "X"))
})

test_that("check_and_fix_chrom_names errors on inconsistent naming", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  ref <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  vcf.df <- data.frame(CHROM = c("chr1", "2"), POS = 1:2)
  expect_error(check_and_fix_chrom_names(vcf.df, ref), "not consistent")
})
