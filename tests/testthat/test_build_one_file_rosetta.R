hepg2_vcf <- function() {
  extdata(
    "annotated-ID-GRCh37",
    "HepG2_Duo_100pM_2mth_cl2.annotated.indel.vcf.gz"
  )
}

# Write a tiny annotated indel VCF to a temp file for targeted tests.
write_mini_vcf <- function(rows) {
  cols <- c(
    "Koh_476", "Koh_89", "COSMIC_83", "long_visual", "ins_or_del_seq",
    "U_seq", "U_seq_count_in_indel_seq", "R", "mh", "unit", "unit_length",
    "internal_rep", "internal_reps", "spacer", "spacer_length", "prime3_rep",
    "prime3_reps", "original_reps"
  )
  df <- as.data.frame(rows, stringsAsFactors = FALSE)
  for (m in setdiff(cols, names(df))) {
    df[[m]] <- NA_character_
  }
  path <- tempfile(fileext = ".annotated.indel.vcf")
  utils::write.table(
    df[, cols],
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  path
}

test_that(".shorten_long_visual keeps flank_5 5' and flank_3 3' bases", {
  x <- "ABCDEFGHIJ <T>[TT] KLMNOPQRSTUVWXYZ"
  expect_equal(
    .shorten_long_visual(x, flank_5 = 3, flank_3 = 4),
    "HIJ <T>[TT] KLMN"
  )
  # Short flanks are kept whole; NA and malformed strings pass through.
  expect_equal(.shorten_long_visual("AB <T> CD", 5, 20), "AB <T> CD")
  expect_equal(.shorten_long_visual(c(NA, "no spaces"), 5, 20),
               c(NA, "no spaces"))
})

test_that(".is_uncapped_del_repeat flags only single-base C/T R10+ classes", {
  expect_true(.is_uncapped_del_repeat("A[Del(C):R12]T"))
  expect_true(.is_uncapped_del_repeat("T[Ins(T):R10]G"))
  expect_false(.is_uncapped_del_repeat("A[Del(C):R9]T"))
  expect_false(.is_uncapped_del_repeat("A[Del(C):R1]T"))
  expect_false(.is_uncapped_del_repeat("Del2:U2:R(5,9)"))
})

test_that(".revcomp_single_tc_long_visual flips only strand-flipped C/T rows", {
  lv <- c(
    "CAGCACTTTG <G>[GG] CGGATCACC",
    "GAGGCTGAGG <G> TGAGTCC",
    "AAA <T> GGG",
    "AAA <TC>[TCTC] GGG",
    NA
  )
  k <- c("G[Del(C):R3]A", "A[Del(C):R1]G", "C[Del(T):R1]A", "Del2:U2:R3", "x")
  out <- .revcomp_single_tc_long_visual(lv, k)
  expect_equal(out[1], "GGTGATCCG <C>[CC] CAAAGTGCTG")
  expect_equal(out[2], "GGACTCA <C> CCTCAGCCTC")
  expect_equal(out[3:5], lv[3:5])
})

test_that("build_one_file_rosetta builds a table from a cell-line VCF", {
  d <- suppressMessages(build_one_file_rosetta(hepg2_vcf()))
  expect_s3_class(d, "data.table")
  expect_equal(
    names(d),
    c("Koh_476", "Koh_89", "COSMIC_83", "n_indels", "example_n", "long_visual")
  )
  expect_gt(nrow(d), 0)

  # Every block has at least one example in context, thanks to the reverse
  # complement of strand-flipped single-base C/T rows.
  expect_false(any(is.na(d$long_visual)))
  expect_false(any(grepl(" ", d$long_visual, fixed = TRUE)))

  # n_indels sums to the number of indels in the VCF.
  raw <- utils::read.delim(gzfile(hepg2_vcf()), check.names = FALSE)
  per_pair <- unique(d[, .(Koh_476, Koh_89, n_indels)])
  expect_equal(sum(per_pair$n_indels), nrow(raw))

  # example_n runs 1..k within each block; blocks are sorted by count.
  expect_true(all(d[, identical(example_n, seq_len(.N)), by = .(Koh_476, Koh_89)]$V1))
  blocks <- d[, .(n = n_indels[1]), by = .(Koh_476, Koh_89)]
  expect_false(is.unsorted(rev(blocks$n)))

  # Labels use open intervals, and 476-types with several 89-types get "*".
  expect_false(any(grepl("R\\(5,9\\)", d$Koh_476)))
  multi <- d[, data.table::uniqueN(Koh_89), by = Koh_476][V1 > 1, Koh_476]
  expect_true(all(endsWith(multi, "*")))
})

test_that("sort_by_count = FALSE orders by the ID476 catalog row order", {
  d <- suppressMessages(build_one_file_rosetta(hepg2_vcf(), sort_by_count = FALSE))
  ord <- catalog_row_order()$ID476
  pos <- match(sub("\\*$", "", d$Koh_476), ord)
  expect_false(any(is.na(pos)))
  expect_false(is.unsorted(pos))
})

test_that("show_details, one_singletc, and n_examples control the columns and rows", {
  d <- suppressMessages(build_one_file_rosetta(
    hepg2_vcf(), show_details = TRUE, one_singletc = TRUE, n_examples = 3
  ))
  expect_true(all(c("ins_or_del_seq", "U_seq", "original_reps") %in% names(d)))
  single <- .is_single_tc_class(sub("\\*$", "", d$Koh_476))
  expect_true(all(d[single, .N, by = .(Koh_476, Koh_89)]$N == 1L))
  expect_true(all(d[!single, .N, by = .(Koh_476, Koh_89)]$N <= 3L))
})

test_that("flank_5 and flank_3 set the flank lengths shown", {
  d <- suppressMessages(build_one_file_rosetta(hepg2_vcf(), flank_5 = 2, flank_3 = 3))
  # Single-base deletions with no remaining repeat: "NN<C>NNN".
  simple <- d[grepl("^[ACGT]{2}<[CT]>[ACGT]{3}$", long_visual)]
  expect_gt(nrow(simple), 0)
  pre <- sub("<.*$", "", d$long_visual)
  post <- sub("^.*[]>}]", "", d$long_visual)
  expect_true(all(nchar(pre) <= 2))
  expect_true(all(nchar(post) <= 3))
})

test_that("cap_9 drops R10+ single-base rows and rewrites R9 as R(9,)", {
  rows <- data.frame(
    Koh_476 = c("A[Del(C):R12]T", "A[Del(C):R9]T", "A[Del(C):R9]T",
                "Del2:U2:R2"),
    Koh_89 = c("A[Del(C):R(9,)]T", "A[Del(C):R(9,)]T", "A[Del(C):R(9,)]T",
               "Del(2,8):U(1,2):R(2,4)"),
    COSMIC_83 = c("DEL:C:1:5+", "DEL:C:1:5+", "DEL:C:1:5+", "DEL:repeats:2:1"),
    long_visual = c(
      "AAAAA <C>[CCCCCCCCCCC] TTTTT",
      "AAAAA <C>[CCCCCCCC] TTTTT",
      "GGGGG <C>[CCCCCCCC] AAAAA",
      "AAAAA <TC>[TC] TTTTT"
    ),
    ins_or_del_seq = c("C", "C", "C", "TC"),
    stringsAsFactors = FALSE
  )
  path <- write_mini_vcf(rows)

  d <- suppressMessages(build_one_file_rosetta(path))
  expect_equal(sort(unique(d$Koh_476)), c("A[Del(C):R(9,)]T", "Del2:U2:R2"))
  expect_equal(unique(d[Koh_476 == "A[Del(C):R(9,)]T", n_indels]), 2L)
  expect_equal(nrow(d[Koh_476 == "A[Del(C):R(9,)]T"]), 2L)

  d0 <- suppressMessages(suppressWarnings(
    build_one_file_rosetta(path, cap_9 = FALSE)
  ))
  expect_true("A[Del(C):R12]T" %in% d0$Koh_476)
  expect_equal(unique(d0[Koh_476 == "A[Del(C):R9]T", n_indels]), 2L)
})

test_that("build_one_file_rosetta writes an xlsx with merged blocks", {
  skip_if_not_installed("openxlsx2")
  out <- tempfile(fileext = ".xlsx")
  d <- suppressMessages(build_one_file_rosetta(hepg2_vcf(), out_path = out))
  expect_true(file.exists(out))

  x <- openxlsx2::read_xlsx(out)
  expect_equal(
    names(x),
    c("476-type", "89-type", "83-type", "N indels", "Example n",
      "Indel in context")
  )
  expect_equal(nrow(x), nrow(d))
  # Values are present in every row (merging does not blank them) ...
  expect_equal(x[["N indels"]], d$n_indels)
  expect_equal(x[["Example n"]], d$example_n)
  # ... and each multi-row block is merged in the four block columns.
  wb <- openxlsx2::wb_load(out)
  merges <- wb$worksheets[[1]]$mergeCells
  n_multi <- d[, .N, by = .(Koh_476, Koh_89)][N > 1, .N]
  expect_equal(length(merges), 4L * n_multi)
  # Rich-text cells drop [] and {} but keep <> around the indel.
  expect_equal(
    x[["Indel in context"]],
    gsub("[][{}]", "", d$long_visual)
  )
})
