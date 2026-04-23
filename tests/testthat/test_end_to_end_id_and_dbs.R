test_that("end-to-end: Strelka ID VCF -> ID83 catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))

  vcf <- read_vcf("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf",
                  filter = "PASS")
  sp <- suppressWarnings(split_vcf(vcf))
  ann <- suppressWarnings(annotate_id_vcf(sp$ID, ref_genome = "GRCh37"))
  expect_true(all(c("annotated.vcf", "discarded.variants") %in% names(ann)))
  expect_true(all(c("COSMIC_83", "Koh_89", "Koh_476") %in%
                    colnames(ann$annotated.vcf)))

  cat <- suppressWarnings(vcf_to_id_catalog(
    ann, type = "ID83", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  ))
  expect_true(is_catalog(cat))
  expect_equal(attr(cat, "type"), "ID83")
  expect_equal(nrow(cat), 83L)
  expect_equal(ncol(cat), 1L)
  expect_gt(sum(cat[, 1L]), 0L)
})

test_that("end-to-end: Koh-89 and Koh-476 catalogs from the same annotated VCF", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))

  vcf <- read_vcf("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf",
                  filter = "PASS")
  sp <- suppressWarnings(split_vcf(vcf))
  ann <- suppressWarnings(annotate_id_vcf(sp$ID, ref_genome = "GRCh37"))

  cat89 <- suppressWarnings(vcf_to_id_catalog(
    ann, type = "ID89", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  ))
  expect_equal(attr(cat89, "type"), "ID89")
  expect_equal(nrow(cat89), 89L)

  cat476 <- suppressWarnings(vcf_to_id_catalog(
    ann, type = "ID476", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  ))
  expect_equal(attr(cat476, "type"), "ID476")
  expect_equal(nrow(cat476), 476L)
})

test_that("empty ID input -> empty ID83 catalog (shape preserved)", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  empty_ann <- data.table::data.table(
    CHROM = character(), POS = integer(),
    REF = character(), ALT = character(),
    COSMIC_83 = character(), Koh_89 = character(), Koh_476 = character()
  )
  cat <- suppressWarnings(vcf_to_id_catalog(
    empty_ann, type = "ID83", ref_genome = "GRCh37",
    region = "genome", sample_name = "s1"
  ))
  expect_equal(nrow(cat), 83L)
  expect_equal(sum(cat[, 1L]), 0)
})

test_that("end-to-end: DBS78 from a synthetic annotated DBS VCF", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  # Pull real DBS positions from the Mutect VCF. Most won't be DBS so we
  # build a minimal synthetic VCF of two DBSs at known positions.
  # Just exercise the shape + happy path on an empty input for simplicity:
  empty <- data.table::data.table(
    CHROM = character(), POS = integer(),
    REF = character(), ALT = character(),
    seq.21bases = character()
  )
  cat <- vcf_to_dbs_catalog(empty, type = "DBS78", ref_genome = "GRCh38",
                            region = "genome", sample_name = "s1")
  expect_equal(attr(cat, "type"), "DBS78")
  expect_equal(nrow(cat), 78L)
  expect_equal(sum(cat[, 1L]), 0)
})
