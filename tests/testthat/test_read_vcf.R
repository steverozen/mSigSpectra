test_that("read_vcf (fread) reads a Mutect VCF and drops ## header lines", {
  out <- read_vcf("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf")
  expect_s3_class(out, "data.table")
  required <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  expect_true(all(required %in% colnames(out)))
  expect_type(out$CHROM, "character")
  expect_type(out$POS, "integer")
  expect_gt(nrow(out), 0L)
})

test_that("read_vcf default filter=TRUE keeps PASS / . / empty rows", {
  unfiltered <- read_vcf("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
                         filter = FALSE)
  filtered <- read_vcf("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
                       filter = TRUE)
  expect_lte(nrow(filtered), nrow(unfiltered))
  expect_true(all(filtered$FILTER %in% c("PASS", ".", "")))
})

test_that("read_vcf accepts a custom filter character vector", {
  passing <- read_vcf("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
                      filter = "PASS")
  expect_true(all(passing$FILTER == "PASS"))
})

test_that("read_vcf reads a Strelka SBS VCF (FILTER-agnostic)", {
  out <- read_vcf("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf")
  expect_gt(nrow(out), 0L)
  expect_true(all(nchar(out$REF) == 1L & nchar(out$ALT) == 1L))
})

test_that("read_vcf reads a Strelka ID VCF", {
  out <- read_vcf("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")
  expect_gt(nrow(out), 0L)
  expect_true(any(nchar(out$REF) != nchar(out$ALT)))
})

test_that("read_vcf backend='vcfppR' errors informatively when vcfppR missing", {
  skip_if(system.file(package = "vcfppR") != "")
  expect_error(
    read_vcf("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
             backend = "vcfppR"),
    "vcfppR"
  )
})

test_that("read_vcf errors on a non-VCF input file", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c("this is not a VCF"), tmp)
  on.exit(unlink(tmp))
  expect_error(read_vcf(tmp), "does not appear to be a VCF")
})

test_that("read_vcfs returns a named list of data.tables", {
  files <- c(
    "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
    "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
  )
  out <- read_vcfs(files)
  expect_length(out, 2L)
  expect_named(out, c("Mutect.GRCh37.s1", "Strelka.SBS.GRCh37.s1"))
  expect_s3_class(out[[1]], "data.table")
  expect_s3_class(out[[2]], "data.table")
})
