test_that("is_grch37 / is_grch38 / is_grcm38 recognize the three shipped genomes", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  skip_if("" == system.file(package = "BSgenome.Mmusculus.UCSC.mm10"))
  expect_true(is_grch37("GRCh37"))
  expect_true(is_grch37("hg19"))
  expect_false(is_grch37("hg38"))
  expect_true(is_grch38("hg38"))
  expect_true(is_grcm38("mm10"))
  expect_false(is_grcm38("GRCh38"))
  expect_false(is_grch37(NULL))
  expect_false(is_grch38(NULL))
})

test_that("infer_trans_ranges returns shipped tables by ref_genome", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  skip_if("" == system.file(package = "BSgenome.Mmusculus.UCSC.mm10"))
  tr <- infer_trans_ranges("GRCh37")
  expect_s3_class(tr, "data.table")
  expect_equal(colnames(tr),
               c("chrom", "start", "end", "strand",
                 "Ensembl.gene.ID", "gene.symbol"))
  expect_identical(tr, mSigSpectra::trans.ranges.GRCh37)
  expect_identical(infer_trans_ranges("hg38"), mSigSpectra::trans.ranges.GRCh38)
  expect_identical(infer_trans_ranges("mm10"), mSigSpectra::trans.ranges.GRCm38)
})

test_that("infer_trans_ranges respects user-supplied trans_ranges", {
  user_tr <- data.table::data.table(chrom = "1", start = 1L, end = 100L,
                                    strand = "+",
                                    Ensembl.gene.ID = "X",
                                    gene.symbol = "Y")
  expect_identical(infer_trans_ranges("GRCh38", trans_ranges = user_tr),
                   user_tr)
})

test_that("add_transcript_strand returns input unchanged for an empty df", {
  vcf <- data.frame(CHROM = character(), POS = integer(), ALT = character())
  out <- add_transcript_strand(vcf, ref_genome = "GRCh38")
  expect_equal(nrow(out), 0L)
})

test_that("add_transcript_strand annotates variants within a known transcript", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  tr <- mSigSpectra::trans.ranges.GRCh38
  first_gene <- tr[1, ]
  vcf <- data.frame(
    CHROM = as.character(first_gene$chrom),
    POS   = as.integer(first_gene$start + 10L),
    REF   = "A", ALT = "G"
  )
  out <- add_transcript_strand(vcf, ref_genome = "GRCh38")
  expect_true("trans.strand" %in% colnames(out))
  expect_equal(as.character(out$trans.strand), as.character(first_gene$strand))
  expect_equal(as.character(out$trans.gene.symbol),
               as.character(first_gene$gene.symbol))
})

test_that("add_transcript_strand leaves unmatched rows with NA strand", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  # A position chosen far from any transcript (coord 1 on chr1)
  vcf <- data.frame(CHROM = "1", POS = 1L, REF = "N", ALT = "A")
  out <- add_transcript_strand(vcf, ref_genome = "GRCh38")
  expect_true(is.na(out$trans.strand))
})
