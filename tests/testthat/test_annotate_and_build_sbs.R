test_that("end-to-end: Strelka SBS VCF -> SBS96 catalog (caller-agnostic)", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))

  # ICAMS's Strelka-specific path reports 770 SBS (its adjacency-merge step
  # folds ~28 adjacent SBS pairs into ~14 DBS rows). mSigSpectra drops that
  # caller-specific merge in favor of a single caller-agnostic path, so the
  # SBS count here is 798 = 770 + 28. That is the correct mSigSpectra
  # output; it is not a bug vs the ICAMS oracle.
  vcf <- read_vcf("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf",
                  filter = "PASS")
  sp <- suppressWarnings(split_vcf(vcf))
  ann <- annotate_vcf(sp$SBS, ref_genome = "GRCh37", variant_type = "SBS")
  cat96 <- vcf_to_catalog(
    ann, type = "SBS96", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  )

  expect_true(is_catalog(cat96))
  expect_equal(attr(cat96, "type"), "SBS96")
  expect_equal(nrow(cat96), 96L)
  expect_equal(ncol(cat96), 1L)
  expect_equal(sum(cat96[, 1L]), 798)   # 770 ICAMS + 28 un-merged adjacencies
})

test_that("end-to-end: Mutect VCF (Mutect.GRCh37.s1) -> SBS96 = 1738 mutations", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))

  vcf <- read_vcf("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
                  filter = "PASS")
  sp <- suppressWarnings(split_vcf(vcf))
  ann <- annotate_vcf(sp$SBS, ref_genome = "GRCh37", variant_type = "SBS")
  cat96 <- vcf_to_catalog(
    ann, type = "SBS96", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  )
  expect_equal(sum(cat96[, 1L]), 1738)
})

test_that("empty SBS VCF -> empty SBS96 catalog with correct shape", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  vcf <- data.frame(
    CHROM = character(), POS = integer(),
    REF = character(), ALT = character(),
    stringsAsFactors = FALSE
  )
  ann <- annotate_vcf(vcf, ref_genome = "GRCh37", variant_type = "SBS")
  # Empty input: vcf_to_catalog should handle by returning a zero-count matrix
  # but it needs an input VCF with the seq.<N>bases column, so construct
  # an empty properly-shaped data.table manually:
  empty_ann <- data.table::data.table(
    CHROM = character(), POS = integer(),
    REF = character(), ALT = character(),
    seq.21bases = character()
  )
  cat96 <- vcf_to_catalog(
    empty_ann, type = "SBS96", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  )
  expect_equal(sum(cat96[, 1L]), 0)
  expect_equal(nrow(cat96), 96L)
})

test_that("pyr_penta canonicalizes A/G-centered pentanucleotides", {
  # ATGCT> T should NOT be reverse-complemented (center C is a pyrimidine)
  # but actually center is G -> flip. Expected: revc("ATGCT") + revc("T") =
  # "AGCAT" + "A" = "AGCATA"
  expect_equal(pyr_penta("ATGCTT"), "AGCATA")
  # ATCGT -> center is C (pyrimidine), unchanged; input "ATCGTA"
  expect_equal(pyr_penta("ATCGTA"), "ATCGTA")
  # Vector
  expect_equal(
    pyr_penta(c("ATGCTT", "ATCGTA", "AAAAAA")),
    c("AGCATA", "ATCGTA", "TTTTTT")
  )
})

test_that("revc_sbs96 reverse-complements a 4-char SBS96 class string", {
  # For "AATC": context="AAT" -> revc="ATT"; target="C" -> revc="G"; "ATT"+"G"
  expect_equal(revc_sbs96("AATC"), "ATTG")
  # For "CCGA": context="CCG" -> revc="CGG"; target="A" -> revc="T"; "CGG"+"T"
  expect_equal(revc_sbs96(c("AATC", "CCGA")), c("ATTG", "CGGT"))
})

test_that("revc returns empty character() for empty input", {
  expect_equal(revc(character(0)), character(0))
})

test_that("vcfs_to_catalogs batches two SBS VCFs into a 2-column SBS96 catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- c(
    "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf",
    "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
  )
  out <- vcfs_to_catalogs(
    files, types = "SBS96",
    ref_genome = "GRCh37", region = "genome"
  )  # Uses read_vcf default filter=TRUE (PASS + . + empty)
  # Not bit-identical to the ICAMS caller-specific filter; hence sums below
  # are the mSigSpectra-agnostic result, not ICAMS's 770.
  expect_named(out, "SBS96")
  expect_equal(ncol(out$SBS96), 2L)
  expect_equal(attr(out$SBS96, "type"), "SBS96")
  expect_equal(colSums(out$SBS96), c(Strelka.SBS.GRCh37.s1 = 798,
                                     Mutect.GRCh37.s1 = 1738))
})
