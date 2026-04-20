test_that("check_and_remove_discarded_variants: empty VCF returns empty df", {
  vcf <- data.frame(CHROM = character(), POS = integer(),
                    REF = character(), ALT = character())
  out <- check_and_remove_discarded_variants(vcf)
  expect_equal(nrow(out$df), 0L)
  expect_null(out$discarded.variants)
})

test_that("rows with REF == ALT are discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1", "2"),
    POS   = c(100L, 200L, 300L),
    REF   = c("A", "C", "G"),
    ALT   = c("A", "T", "G"),   # row 1 and 3 have REF==ALT
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 1L)
  expect_equal(nrow(out$discarded.variants), 2L)
  expect_true(all(out$discarded.variants$discarded.reason ==
                    "Variant with same REF and ALT"))
})

test_that("multi-ALT rows are discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1"), POS = c(100L, 200L),
    REF   = c("A", "C"), ALT = c("G,T", "T"),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 1L)
  expect_true("Variant with multiple alternative alleles"
              %in% out$discarded.variants$discarded.reason)
})

test_that("exact duplicate rows collapse to one", {
  vcf <- data.frame(
    CHROM = c("1", "1", "1"),
    POS   = c(100L, 100L, 200L),
    REF   = c("A", "A", "C"),
    ALT   = c("G", "G", "T"),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 2L)
  expect_equal(nrow(out$discarded.variants), 1L)
})

test_that("duplicate CHROM/POS/REF with different ALT are all discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1", "2"),
    POS   = c(100L, 100L, 200L),
    REF   = c("A", "A", "C"),
    ALT   = c("G", "T", "G"),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 1L)    # only the chrom 2 variant survives
  expect_equal(nrow(out$discarded.variants), 2L)
})

test_that("3+ bp substitutions are discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1"), POS = c(100L, 200L),
    REF   = c("A", "ACG"), ALT = c("T", "TGA"),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 1L)
  expect_true("Variant involves three or more nucleotides"
              %in% out$discarded.variants$discarded.reason)
})

test_that("complex indels (REF[1] != ALT[1]) are discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1"), POS = c(100L, 200L),
    REF   = c("A", "ATG"), ALT = c("ATG", "CGG"),  # A -> ATG is clean insertion
    stringsAsFactors = FALSE                        # ATG -> CGG is complex
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  # ATG -> CGG: REF[1]='A' != ALT[1]='C' AND lengths equal -> caught by "3+ bp subs"
  # A -> ATG: REF[1]='A' == ALT[1]='A' -> passes as clean insertion
  expect_equal(nrow(out$df), 1L)
  expect_equal(out$df$REF, "A")
})

test_that("ambiguous REF bases are discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1"), POS = c(100L, 200L),
    REF   = c("N", "A"), ALT = c("A", "G"),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 1L)
  expect_true("Variant has ambiguous REF base"
              %in% out$discarded.variants$discarded.reason)
})

test_that("wrong DBS (shared base) rows are discarded", {
  vcf <- data.frame(
    CHROM = c("1", "1"), POS = c(100L, 200L),
    REF   = c("TA", "AT"), ALT = c("TT", "GC"),  # 'TA>TT' shares position 1
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(check_and_remove_discarded_variants(vcf))
  expect_equal(nrow(out$df), 1L)
  expect_true("Wrong DBS variant"
              %in% out$discarded.variants$discarded.reason)
})
