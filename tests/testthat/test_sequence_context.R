test_that("add_seq_context attaches a seq.<N>bases column of the expected width", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  vcf <- data.frame(CHROM = c("1", "1"), POS = c(1000000L, 2000000L))
  out <- add_seq_context(vcf, ref_genome = "GRCh37", seq_context_width = 10)
  expect_true("seq.21bases" %in% colnames(out))
  expect_equal(nchar(out[["seq.21bases"]]), c(21L, 21L))
  expect_match(out[["seq.21bases"]], "^[ACGTN]{21}$")
})

test_that("add_seq_context respects a custom seq_context_width", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  vcf <- data.frame(CHROM = "1", POS = 5000000L)
  out <- add_seq_context(vcf, ref_genome = "GRCh38", seq_context_width = 5)
  expect_true("seq.11bases" %in% colnames(out))
  expect_equal(nchar(out[["seq.11bases"]]), 11L)
})

test_that("add_seq_context returns input unchanged for an empty df", {
  vcf <- data.frame(CHROM = character(), POS = integer())
  out <- add_seq_context(vcf, ref_genome = "GRCh38")
  expect_equal(nrow(out), 0L)
})

test_that("add_seq_context harmonizes chr-prefix mismatches", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  vcf <- data.frame(CHROM = "1", POS = 5000000L)  # no "chr" prefix; hg38 uses "chr"
  out <- add_seq_context(vcf, ref_genome = "GRCh38", seq_context_width = 3)
  expect_equal(nchar(out[["seq.7bases"]]), 7L)
})

test_that("normalize_genome_arg recognizes common aliases", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  g1 <- normalize_genome_arg("hg38")
  g2 <- normalize_genome_arg("GRCh38")
  g3 <- normalize_genome_arg("BSgenome.Hsapiens.UCSC.hg38")
  expect_identical(BSgenome::organism(g1), BSgenome::organism(g2))
  expect_identical(BSgenome::organism(g1), BSgenome::organism(g3))
})

test_that("normalize_genome_arg errors on unknown identifiers", {
  expect_error(normalize_genome_arg("invalid_name"), "Unrecognized")
  expect_error(normalize_genome_arg(NULL), "non-NULL")
})

test_that("stop_if_region_illegal accepts valid regions and errors on others", {
  expect_silent(stop_if_region_illegal("genome"))
  expect_silent(stop_if_region_illegal("exome"))
  expect_silent(stop_if_region_illegal("transcript"))
  expect_silent(stop_if_region_illegal("unknown"))
  expect_error(stop_if_region_illegal("wholegenome"), "Unrecognized")
})

test_that("stop_if_transcribed_region_illegal rejects genome region", {
  expect_error(stop_if_transcribed_region_illegal("genome"),
               "transcript, exome, or unknown")
  expect_silent(stop_if_transcribed_region_illegal("transcript"))
})
