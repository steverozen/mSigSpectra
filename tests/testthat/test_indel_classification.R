test_that("justify_indel moves a deletion left as far as possible", {
  r <- justify_indel("CAAAG", "CAAG", pos = 2, expected_delta = "A")
  expect_equal(r$leftmost_pos, 2L)
  expect_null(r$error)

  # Same deletion declared at a different position should still normalize
  r2 <- justify_indel("CAAAG", "CAAG", pos = 3, expected_delta = "A")
  expect_equal(r2$leftmost_pos, 2L)
  expect_null(r2$error)
})

test_that("justify_indel handles a 2bp deletion in a repeat", {
  r <- justify_indel("CACAG", "CAG", pos = 3, expected_delta = "CA")
  expect_equal(r$leftmost_pos, 2L)
  expect_equal(r$del_str, "AC")  # leftmost delete is "AC" starting at pos 2
  expect_true(r$edge_warning)  # ran off left edge
})

test_that("justify_indel reports error on mismatched expected_delta", {
  r <- justify_indel("CAAAG", "CAAG", pos = 2, expected_delta = "G")
  expect_false(is.null(r$error))
  expect_match(r$error, "expected_delta")
})

test_that("categorize_1_justified_indel: 1bp T/A deletion in homopolymer -> DEL:T:1:R", {
  # GGAAAGG: 3 A's; deleting one A at pos 3 leaves 2 A's -> R = 2 (counted
  # as 2 repeat units after the deletion)
  r <- categorize_1_justified_indel("GGAAAGG", "d", ins_or_del_seq = "A", pos = 3)
  expect_equal(r$COSMIC_83, "DEL:T:1:2")
  expect_equal(r$ins_or_del, "d")
  # Note: categorize_1_justified_indel internally reverse-complements A to T
  # so ins_or_del_seq in the return value is "T"
  expect_equal(r$ins_or_del_seq, "T")
})

test_that("categorize_1_justified_indel: 1bp C/G deletion -> DEL:C:*:*", {
  r <- categorize_1_justified_indel("AAACCCAA", "d", ins_or_del_seq = "C", pos = 4)
  expect_true(startsWith(r$COSMIC_83, "DEL:C:1:"))
})

test_that("categorize_1_justified_indel: 1bp insertion", {
  r <- categorize_1_justified_indel("AATTAA", "i", ins_or_del_seq = "T", pos = 4)
  expect_equal(r$ins_or_del, "i")
  expect_true(startsWith(r$COSMIC_83, "INS:T:1:"))
})

test_that("lcprefix_fast returns length of matching prefix", {
  expect_equal(lcprefix_fast("ACGT", "ACCA"), 2L)
  expect_equal(lcprefix_fast("ACGT", "ACGT"), 4L)
  expect_equal(lcprefix_fast("ACGT", "TCGA"), 0L)
  expect_equal(lcprefix_fast("", "ACGT"), 0L)
  expect_equal(lcprefix_fast("AAAA", "AA"), 2L)
})

test_that("seg_simple returns the expected structure for an AT-repeat deletion", {
  r <- seg_simple("d", "ATAT", "ATATGG")
  expect_type(r, "list")
  expect_named(r, c("unit", "unit_length", "internal_rep", "internal_reps",
                    "spacer", "spacer_length", "prime3_rep", "prime3_reps",
                    "original_reps"))
  expect_equal(r$unit, "AT")
})

test_that("justify_id_vcf end-to-end on a real Strelka ID VCF", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  vcf <- read_vcf("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf",
                  filter = "PASS")
  # Keep only the first 20 rows to keep this test fast.
  vcf <- vcf[seq_len(min(20L, nrow(vcf))), ]
  out <- suppressWarnings(justify_id_vcf(vcf, ref.genome = "GRCh37",
                                         explain_indels = 0))
  expect_type(out, "list")
  expect_true("annotated.vcf" %in% names(out))
  expect_true("seq.context" %in% colnames(out$annotated.vcf))
  expect_true("pos_shift" %in% colnames(out$annotated.vcf))
})

test_that("annot_vcf_to_83_catalog tallies COSMIC_83 strings into a 1-col ID83 matrix", {
  # Build a synthetic VCF with COSMIC_83 column already populated.
  vcf <- data.frame(
    CHROM = rep("1", 4),
    POS   = c(1L, 2L, 3L, 4L),
    REF   = c("A", "C", "G", "T"),
    ALT   = c("A", "C", "G", "T"),
    COSMIC_83 = c("DEL:T:1:2", "DEL:T:1:2", "INS:C:1:0", "DEL:C:1:1"),
    stringsAsFactors = FALSE
  )
  out <- annot_vcf_to_83_catalog(vcf, sample_id = "s1")
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 83L)
  expect_equal(ncol(out), 1L)
  expect_equal(sum(out[, 1L]), 4L)
  expect_equal(out["DEL:T:1:2", "s1"], 2)
  expect_equal(out["INS:C:1:0", "s1"], 1)
  expect_equal(out["DEL:C:1:1", "s1"], 1)
})
