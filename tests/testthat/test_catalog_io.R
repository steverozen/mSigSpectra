make_sbs96_cat <- function() {
  rns <- mSigSpectra::catalog.row.order$SBS96
  m <- matrix(seq_len(96L), nrow = 96L, ncol = 1L,
              dimnames = list(rns, "sample1"))
  as_catalog(m, ref_genome = NULL, region = "genome")
}

test_that("write_catalog + read_catalog round-trip SBS96 losslessly (ICAMS)", {
  cat_in <- make_sbs96_cat()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write_catalog(cat_in, tmp, format = "ICAMS")
  cat_out <- read_catalog(tmp, region = "genome")
  expect_equal(attr(cat_out, "type"), "SBS96")
  expect_equal(as.numeric(cat_out[, 1L]), as.numeric(cat_in[, 1L]))
  expect_equal(rownames(cat_out), rownames(cat_in))
})

test_that("round-trip DBS78", {
  rns <- mSigSpectra::catalog.row.order$DBS78
  m <- matrix(seq_len(78L), nrow = 78L, ncol = 1L,
              dimnames = list(rns, "sample1"))
  cat_in <- as_catalog(m, ref_genome = NULL, region = "genome")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write_catalog(cat_in, tmp)
  cat_out <- read_catalog(tmp, region = "genome")
  expect_equal(attr(cat_out, "type"), "DBS78")
  expect_equal(as.numeric(cat_out[, 1L]), as.numeric(cat_in[, 1L]))
  expect_equal(rownames(cat_out), rownames(cat_in))
})

test_that("round-trip ID83", {
  rns <- mSigSpectra::catalog.row.order$ID
  m <- matrix(seq_len(83L), nrow = 83L, ncol = 1L,
              dimnames = list(rns, "sample1"))
  cat_in <- as_catalog(m, ref_genome = NULL, region = "genome")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write_catalog(cat_in, tmp)
  cat_out <- read_catalog(tmp, region = "genome")
  expect_equal(attr(cat_out, "type"), "ID83")
  expect_equal(as.numeric(cat_out[, 1L]), as.numeric(cat_in[, 1L]))
  expect_equal(rownames(cat_out), rownames(cat_in))
})

test_that("round-trip SBS1536", {
  rns <- mSigSpectra::catalog.row.order$SBS1536
  m <- matrix(seq_len(1536L), nrow = 1536L, ncol = 1L,
              dimnames = list(rns, "sample1"))
  cat_in <- as_catalog(m, ref_genome = NULL, region = "genome")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write_catalog(cat_in, tmp)
  cat_out <- read_catalog(tmp, region = "genome")
  expect_equal(attr(cat_out, "type"), "SBS1536")
  expect_equal(sum(cat_out[, 1L]), sum(cat_in[, 1L]))
})

test_that("detect_catalog_format identifies SigProfiler vs ICAMS", {
  sigpro_dt <- data.table::data.table(
    mut = c("A[C>A]A", "A[C>A]C"),
    s1 = c(1, 2)
  )
  expect_equal(detect_catalog_format(sigpro_dt), "SigProfiler")

  icams_dt <- data.table::data.table(
    `Mutation type` = c("C>A", "C>A"),
    Trinucleotide = c("ACA", "ACC"),
    s1 = c(1, 2)
  )
  expect_equal(detect_catalog_format(icams_dt), "ICAMS")
})

test_that("write_catalog errors on SigProfiler / COSMIC formats (not yet implemented)", {
  cat_in <- make_sbs96_cat()
  tmp <- tempfile(fileext = ".csv")
  expect_error(write_catalog(cat_in, tmp, format = "SigProfiler"),
               "not yet implemented")
  expect_error(write_catalog(cat_in, tmp, format = "COSMIC"),
               "not yet implemented")
})

test_that("read_catalog errors on an incompatible catalog file", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c("foo,bar,s1", "x,y,1"), tmp)
  expect_error(read_catalog(tmp), class = "error")
})

test_that("read_catalog (SigProfiler SBS96) parses A[C>A]A labels", {
  # Build a synthetic SigProfiler-style SBS96 file
  labels <- c(sapply(c("A", "C", "G", "T"), function(a) {
    sapply(c("A", "C", "G", "T"), function(b) {
      sapply(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), function(m) {
        paste0(a, "[", m, "]", b)
      })
    })
  }))
  stopifnot(length(labels) == 96L)
  dt <- data.table::data.table(mut = labels, s1 = seq_len(96L))
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  data.table::fwrite(dt, tmp, sep = "\t")
  out <- read_catalog(tmp, region = "genome", format = "SigProfiler")
  expect_equal(attr(out, "type"), "SBS96")
  expect_equal(sum(out[, 1L]), sum(seq_len(96L)))
})
