inst_dir <- system.file(package = "mSigSpectra")

skip_reijns <- function() {
  skip_on_cran()
  skip_if(
    "" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"),
    "BSgenome.Hsapiens.UCSC.hg38 not installed"
  )
  skip_if(
    length(list.files(inst_dir, pattern = "reijns_cell.*\\.vcf$")) == 0L,
    "reijns_cell*.vcf files not present in inst/"
  )
}

source(system.file("make_reijns_rosetta_csv.R", package = "mSigSpectra"),
       local = TRUE)
source(system.file("make_reijns_catalog.R",     package = "mSigSpectra"),
       local = TRUE)

# ---- rosetta CSV snapshots ---------------------------------------------------

test_that("make_reijns_rosetta_csv keep_high_R snapshot", {
  skip_reijns()
  out <- withr::local_tempfile(fileext = ".csv")
  make_reijns_rosetta_csv(
    vcf_dir     = inst_dir,
    ref_genome  = "hg38",
    out_path    = out,
    drop_high_R = FALSE
  )
  expect_snapshot_file(out, "reijns_rosetta_keep_high_R.csv")
})

test_that("make_reijns_rosetta_csv drop_high_R snapshot", {
  skip_reijns()
  out <- withr::local_tempfile(fileext = ".csv")
  make_reijns_rosetta_csv(
    vcf_dir     = inst_dir,
    ref_genome  = "hg38",
    out_path    = out,
    drop_high_R = TRUE
  )
  expect_snapshot_file(out, "reijns_rosetta_drop_high_R.csv")
})

# ---- catalog TSV snapshots ---------------------------------------------------

test_that("make_reijns_catalog keep_high_R snapshot", {
  skip_reijns()
  out_keep <- withr::local_tempfile(fileext = ".tsv")
  out_drop <- withr::local_tempfile(fileext = ".tsv")
  make_reijns_catalog(
    vcf_dir    = inst_dir,
    ref_genome = "hg38",
    out_keep   = out_keep,
    out_drop   = out_drop
  )
  expect_snapshot_file(out_keep, "reijns_catalog_keep_high_R.tsv")
})

test_that("make_reijns_catalog drop_high_R snapshot", {
  skip_reijns()
  out_keep <- withr::local_tempfile(fileext = ".tsv")
  out_drop <- withr::local_tempfile(fileext = ".tsv")
  make_reijns_catalog(
    vcf_dir    = inst_dir,
    ref_genome = "hg38",
    out_keep   = out_keep,
    out_drop   = out_drop
  )
  expect_snapshot_file(out_drop, "reijns_catalog_drop_high_R.tsv")
})
