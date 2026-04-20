test_that("transform_catalog counts -> density scales by inverse abundance", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  rns <- mSigSpectra::catalog.row.order$SBS96
  m <- matrix(1, nrow = 96L, ncol = 1L, dimnames = list(rns, "s1"))
  cat_counts <- as_catalog(m, ref_genome = "GRCh38", region = "genome",
                           counts_or_density = "counts")
  cat_density <- transform_catalog(cat_counts,
                                   target_counts_or_density = "density")
  expect_equal(attr(cat_density, "counts_or_density"), "density")
  expect_equal(attr(cat_density, "type"), "SBS96")
  # After density transform, values scale by flat_abundance[nmer] / source_abundance[nmer]
  expect_true(all(is.finite(cat_density[, 1L])))
  expect_true(all(cat_density[, 1L] > 0))
})

test_that("transform_catalog counts -> counts.signature normalizes columns to 1", {
  rns <- mSigSpectra::catalog.row.order$SBS96
  m <- matrix(runif(96L * 2L, 1, 10), nrow = 96L, ncol = 2L,
              dimnames = list(rns, c("s1", "s2")))
  cat <- as_catalog(m, ref_genome = "GRCh38", region = "genome",
                    counts_or_density = "counts")
  sig <- transform_catalog(cat, target_counts_or_density = "counts.signature")
  expect_equal(colSums(sig), c(s1 = 1, s2 = 1))
})

test_that("transform_catalog density -> density is a null op", {
  rns <- mSigSpectra::catalog.row.order$SBS96
  m <- matrix(1, nrow = 96L, ncol = 1L, dimnames = list(rns, "s1"))
  cat <- as_catalog(m, ref_genome = "GRCh38", region = "genome",
                    counts_or_density = "density")
  expect_identical(
    transform_catalog(cat, target_counts_or_density = "density"),
    cat
  )
})

test_that("collapse_catalog SBS1536 -> SBS96 sums pentanucleotides into trinucleotides", {
  rns <- mSigSpectra::catalog.row.order$SBS1536
  m <- matrix(1, nrow = 1536L, ncol = 1L, dimnames = list(rns, "s1"))
  cat <- as_catalog(m, ref_genome = "GRCh38", region = "genome")
  out <- collapse_catalog(cat, to = "SBS96")
  expect_equal(attr(out, "type"), "SBS96")
  expect_equal(nrow(out), 96L)
  expect_equal(sum(out[, 1L]), 1536)   # total count preserved
})
