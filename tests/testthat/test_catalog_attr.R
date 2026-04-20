make_test_sbs96 <- function() {
  m <- matrix(1, nrow = 96L, ncol = 2L,
              dimnames = list(mSigSpectra::catalog.row.order$SBS96,
                              c("s1", "s2")))
  m
}

test_that("as_catalog infers type from nrow and attaches attributes", {
  cat1 <- as_catalog(make_test_sbs96())
  expect_equal(attr(cat1, "type"), "SBS96")
  expect_equal(attr(cat1, "counts_or_density"), "counts")
  expect_equal(attr(cat1, "region"), "unknown")
  expect_null(attr(cat1, "class"))   # no S3 class
  expect_true(is.matrix(cat1))
})

test_that("as_catalog infers type for every supported nrow", {
  types <- c("SBS96" = 96L, "SBS192" = 192L, "SBS1536" = 1536L,
             "DBS78" = 78L, "DBS136" = 136L, "DBS144" = 144L,
             "ID83" = 83L, "ID166" = 166L)
  for (nm in names(types)) {
    n <- types[[nm]]
    key <- if (nm == "ID83") "ID" else nm
    rns <- mSigSpectra::catalog.row.order[[key]]
    m <- matrix(1, nrow = n, ncol = 1L, dimnames = list(rns, "s1"))
    # SBS192 / DBS144 need region in {exome, transcript, unknown}; use transcript
    if (nm %in% c("SBS192", "DBS144")) {
      cat <- as_catalog(m, region = "transcript")
    } else {
      cat <- as_catalog(m)
    }
    expect_equal(attr(cat, "type"), nm, info = nm)
  }
})

test_that("as_catalog accepts a named numeric vector and sets colname Unknown", {
  v <- setNames(rep(1, 96L), mSigSpectra::catalog.row.order$SBS96)
  cat <- as_catalog(v)
  expect_equal(ncol(cat), 1L)
  expect_equal(colnames(cat), "Unknown")
  expect_equal(attr(cat, "type"), "SBS96")
})

test_that("as_catalog reorders rows to the canonical order", {
  m <- make_test_sbs96()
  m_shuffled <- m[sample(96L), ]
  cat <- as_catalog(m_shuffled)
  expect_equal(rownames(cat), mSigSpectra::catalog.row.order$SBS96)
})

test_that("as_catalog errors on missing rownames unless infer_rownames=TRUE", {
  m <- matrix(1, nrow = 96L, ncol = 1L)
  expect_error(as_catalog(m), "rownames")
  cat <- as_catalog(m, infer_rownames = TRUE)
  expect_equal(rownames(cat), mSigSpectra::catalog.row.order$SBS96)
})

test_that("as_catalog errors on unsupported nrow", {
  m <- matrix(1, nrow = 100L, ncol = 1L,
              dimnames = list(as.character(1:100), "s1"))
  expect_error(as_catalog(m), "Unsupported number of rows")
})

test_that("as_catalog errors on illegal region", {
  m <- make_test_sbs96()
  expect_error(as_catalog(m, region = "wholegenome"), "Unrecognized region")
})

test_that("as_catalog for SBS192/DBS144 promotes region='genome' to 'transcript'", {
  rns <- mSigSpectra::catalog.row.order$SBS192
  m <- matrix(1, nrow = 192L, ncol = 1L, dimnames = list(rns, "s1"))
  cat <- as_catalog(m, region = "genome")
  expect_equal(attr(cat, "region"), "transcript")
})

test_that("is_catalog recognizes a properly constructed catalog", {
  cat <- as_catalog(make_test_sbs96())
  expect_true(is_catalog(cat))
  expect_false(is_catalog(make_test_sbs96()))   # plain matrix, no attrs
  expect_false(is_catalog("not a matrix"))
})

test_that("catalog_attrs returns all five attribute fields", {
  cat <- as_catalog(make_test_sbs96(), region = "genome")
  a <- catalog_attrs(cat)
  expect_named(a, c("type", "counts_or_density", "ref_genome", "region",
                    "abundance"))
  expect_equal(a$type, "SBS96")
  expect_equal(a$region, "genome")
})

test_that("cbind_catalogs concatenates samples and preserves attributes", {
  a <- as_catalog(make_test_sbs96())
  b <- as_catalog(make_test_sbs96())
  out <- cbind_catalogs(list(a, b))
  expect_equal(ncol(out), 4L)
  expect_equal(attr(out, "type"), "SBS96")
  expect_equal(rownames(out), mSigSpectra::catalog.row.order$SBS96)
})

test_that("cbind_catalogs rejects catalogs of different types", {
  a <- as_catalog(make_test_sbs96())
  m192 <- matrix(1, nrow = 192L, ncol = 1L,
                 dimnames = list(mSigSpectra::catalog.row.order$SBS192, "s1"))
  b <- as_catalog(m192, region = "transcript")
  expect_error(cbind_catalogs(list(a, b)), "mismatched")
})

test_that("subset_catalog preserves attributes (except abundance when rows subset)", {
  cat <- as_catalog(make_test_sbs96(), ref_genome = NULL)
  sub <- subset_catalog(cat, cols = 1L)
  expect_equal(ncol(sub), 1L)
  expect_equal(attr(sub, "type"), "SBS96")
  expect_equal(attr(sub, "region"), "unknown")
})
