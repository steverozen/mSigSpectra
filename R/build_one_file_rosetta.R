#' Build a one-file indel "Rosetta stone" Excel table
#'
#' Reads a single annotated indel VCF (as written by [annotate_id_vcf()]) and
#' writes an Excel workbook that links the `Koh_476`, `Koh_89`, and
#' `COSMIC_83` classifications, with example indels in sequence context for
#' each observed (`Koh_476`, `Koh_89`) pair. This is a single-file version of
#' the three-step pipeline that built Supplementary Table S1 of Liu et al.
#' (2026), with the steps combined into one function.
#'
#' @param vcf_path Path to one annotated indel VCF (may be gzipped). Must
#'   contain the columns `Koh_476`, `Koh_89`, `COSMIC_83`, and `long_visual`.
#' @param out_path Path of the `.xlsx` file to write. If `NULL`, no file is
#'   written and only the table is returned.
#' @param flank_5 Number of bases to keep from the 5' flank in `long_visual`.
#' @param flank_3 Number of bases to keep from the 3' flank in `long_visual`.
#' @param cap_9 If `TRUE`, drop single-base C/T indel rows whose repeat count
#'   is 10 or more (no counterpart in the cap-9 ID476/ID89 mapping) and rewrite
#'   `R9` as `R(9,)`.
#' @param show_details If `TRUE`, include the annotation detail columns after
#'   "Indel in context".
#' @param one_singletc If `TRUE`, show 1 example per single-base C/T class
#'   instead of 5.
#' @param n_examples Maximum number of examples per (`Koh_476`, `Koh_89`) pair
#'   for classes other than single-base C/T indels.
#' @param sort_by_count If `TRUE`, order the (`Koh_476`, `Koh_89`) blocks by
#'   descending `n_indels` (the number of indels in the input VCF with that
#'   pair) instead of by the canonical ID476 row order (ties broken by that
#'   order). LibreOffice and Excel
#'   cannot sort a range containing merged cells, so this is the way to get
#'   a count-sorted table.
#'
#' @return Invisibly, a `data.table` with the rows written to the workbook.
#'
#' @export
build_one_file_rosetta <- function(
  vcf_path,
  out_path = NULL,
  flank_5 = 5,
  flank_3 = 20,
  cap_9 = TRUE,
  show_details = FALSE,
  one_singletc = FALSE,
  n_examples = 20,
  sort_by_count = TRUE
) {
  if (!is.null(out_path) && !requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("Package 'openxlsx2' is required to write the Excel file")
  }

  # ---- Step 1: read the VCF and keep the classification columns ----
  keep_cols <- c(
    "Koh_476",
    "Koh_89",
    "COSMIC_83",
    "long_visual",
    "ins_or_del_seq",
    "U_seq",
    "U_seq_count_in_indel_seq",
    "R",
    "mh",
    "unit",
    "unit_length",
    "internal_rep",
    "internal_reps",
    "spacer",
    "spacer_length",
    "prime3_rep",
    "prime3_reps",
    "original_reps"
  )
  extra_cols <- keep_cols[-(1:4)]

  message("Reading ", basename(vcf_path))
  df <- utils::read.delim(
    gzfile(vcf_path),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  missing <- setdiff(keep_cols, colnames(df))
  if (length(missing) > 0) {
    warning(vcf_path, " is missing columns: ", paste(missing, collapse = ", "))
    for (m in missing) {
      df[[m]] <- NA_character_
    }
  }
  df$long_visual <- .shorten_long_visual(df$long_visual, flank_5, flank_3)
  full <- data.table::as.data.table(df[, keep_cols, drop = FALSE])

  # ---- Step 2: cap-9 filter ----
  if (cap_9) {
    drop <- .is_uncapped_del_repeat(full$Koh_476)
    message("Dropping ", sum(drop), " of ", nrow(full), " uncapped rows")
    full <- full[!drop]
    # Rewrite "(Ins|Del)(C|T):R9]" as "(Ins|Del)(C|T):R(9,)]" so the Koh_476
    # values match the ID476_ID89 mapping convention for the open-ended bin.
    full$Koh_476 <- sub(
      "((?:Ins|Del)\\([CT]\\)):R9\\]",
      "\\1:R(9,)]",
      full$Koh_476,
      perl = TRUE
    )
  }
  pairs <- unique(full[, c("Koh_476", "Koh_89"), with = FALSE])
  data.table::setorder(pairs, Koh_476, Koh_89)
  message("Found ", nrow(pairs), " (Koh_476, Koh_89) pairs")
  # Number of indels (after the cap-9 filter) with each pair.
  pair_counts <- full[, .(n_indels = .N), by = .(Koh_476, Koh_89)]

  # ---- Step 3: pick examples and build the table ----
  n_examples_singletc <- if (one_singletc) 1L else 5L

  # For single-base T/C indels, restrict examples to ones where the indel base
  # on the reported strand is itself C or T (long_visual middle token like
  # "<C>" or "<T>"). n_examples_singletc such examples are kept per class.
  single_tc <- .is_single_tc_class(full$Koh_476)
  drop <- single_tc & !grepl(" <[CT]>", full$long_visual)
  message(
    "Dropping ",
    sum(drop),
    " rows whose strand-flipped base hides the canonical C/T"
  )
  full <- full[!drop]

  full_unique <- unique(full, by = c("Koh_476", "Koh_89", "long_visual"))
  data.table::setorder(full_unique, Koh_476, Koh_89)
  examples <- full_unique[,
    {
      if (.is_single_tc_class(.BY$Koh_476)) {
        utils::head(.SD, n_examples_singletc)
      } else {
        seqs <- unique(ins_or_del_seq)
        if (length(seqs) >= n_examples) {
          # More (or equal) seqs than budget: one row per seq, capped.
          .SD[match(seqs[seq_len(n_examples)], ins_or_del_seq)]
        } else {
          # Fewer seqs than budget: round-robin to fill up to n_examples.
          # seq_rank    = which seq each row belongs to (position in seqs).
          # within_rank = ordinal of the row within its seq group.
          # Sorting by (within_rank, seq_rank) interleaves seqs so every seq
          # gets one example before any seq gets a second.
          seq_rank <- match(ins_or_del_seq, seqs)
          within_rank <- stats::ave(seq_rank, seq_rank, FUN = seq_along)
          utils::head(.SD[order(within_rank, seq_rank)], n_examples)
        }
      }
    },
    by = .(Koh_476, Koh_89),
    .SDcols = c("COSMIC_83", "long_visual", extra_cols)
  ]

  cosmic_per_pair <- full_unique[,
    .(COSMIC_83 = paste(sort(unique(COSMIC_83)), collapse = "; ")),
    by = .(Koh_476, Koh_89)
  ]

  doc <- merge(
    pairs,
    examples[, c("Koh_476", "Koh_89", "long_visual", extra_cols), with = FALSE],
    by = c("Koh_476", "Koh_89"),
    all.x = TRUE,
    sort = FALSE
  )
  doc <- merge(
    doc,
    cosmic_per_pair,
    by = c("Koh_476", "Koh_89"),
    all.x = TRUE,
    sort = FALSE
  )
  doc <- merge(
    doc,
    pair_counts,
    by = c("Koh_476", "Koh_89"),
    all.x = TRUE,
    sort = FALSE
  )
  doc[, example_n := seq_len(.N), by = .(Koh_476, Koh_89)]

  # Convert to open-interval labels first, because catalog_row_order()$ID476
  # uses them (e.g. "R(5,)" rather than the older "R(5,9)").
  doc[, Koh_476 := change_476_type_ids_to_open_intervals(Koh_476)]
  doc[, Koh_89 := change_89_type_ids_to_open_intervals(Koh_89)]

  # Order by the canonical ID476 row order.
  ord <- catalog_row_order()$ID476
  doc[, .row_ord := match(Koh_476, ord)]
  not_in_order <- unique(doc$Koh_476[is.na(doc$.row_ord)])
  if (length(not_in_order) > 0) {
    warning(
      "Koh_476 values not in catalog_row_order()$ID476: ",
      paste(not_in_order, collapse = ", ")
    )
  }
  if (sort_by_count) {
    data.table::setorder(
      doc,
      -n_indels,
      .row_ord,
      Koh_89,
      example_n,
      na.last = TRUE
    )
  } else {
    data.table::setorder(doc, .row_ord, Koh_89, example_n, na.last = TRUE)
  }
  doc[, .row_ord := NULL]

  # Collapse spaces in every long_visual.
  doc$long_visual <- gsub(" ", "", doc$long_visual, fixed = TRUE)

  # Asterisk Koh_476 values that map to more than one Koh_89.
  multi_koh476 <- doc[, data.table::uniqueN(Koh_89), by = Koh_476][
    V1 > 1,
    Koh_476
  ]
  if (length(multi_koh476) > 0) {
    message(
      "Koh_476 classes mapping to multiple Koh_89 (marked with *): ",
      paste(multi_koh476, collapse = ", ")
    )
    doc[Koh_476 %in% multi_koh476, Koh_476 := paste0(Koh_476, "*")]
  }

  doc <- doc[,
    c(
      "Koh_476",
      "Koh_89",
      "COSMIC_83",
      "n_indels",
      "example_n",
      "long_visual",
      extra_cols
    ),
    with = FALSE
  ]
  if (!show_details) {
    doc <- doc[, seq_len(which(names(doc) == "long_visual")), with = FALSE]
  }

  n_pairs <- data.table::uniqueN(doc, by = c("Koh_476", "Koh_89"))
  counts <- doc[, .N, by = .(Koh_476, Koh_89)]
  n_thin <- counts[N < n_examples, .N]
  message("Pairs written: ", n_pairs)
  message("Pairs with < ", n_examples, " examples: ", n_thin)

  if (!is.null(out_path)) {
    .write_rosetta_xlsx(doc, out_path)
    message("Wrote ", out_path)
  }
  invisible(doc)
}

# long_visual is "<5'-flank> <indel-token> <3'-flank>". Keep the last flank_5
# chars of the 5' flank, the indel token, and the first flank_3 chars of the
# 3' flank.
.shorten_long_visual <- function(x, flank_5 = 5, flank_3 = 20) {
  out <- x
  ok <- !is.na(x)
  parts <- strsplit(x[ok], " ", fixed = TRUE)
  ok_three <- vapply(parts, length, integer(1)) == 3L
  idx <- which(ok)[ok_three]
  good <- parts[ok_three]
  pre <- vapply(good, `[`, character(1), 1)
  mid <- vapply(good, `[`, character(1), 2)
  post <- vapply(good, `[`, character(1), 3)
  out[idx] <- paste(
    substr(pre, pmax(1, nchar(pre) - flank_5 + 1), nchar(pre)),
    mid,
    substr(post, 1, flank_3)
  )
  out
}

# Returns TRUE for Koh_476 strings like "X[Del(C):Rnn]Y", "X[Del(T):Rnn]Y",
# "X[Ins(C):Rnn]Y", or "X[Ins(T):Rnn]Y" where nn is a repeat count of 10+.
.is_uncapped_del_repeat <- function(x) {
  grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R1[0-9]+\\][ACGT]$", x)
}

.is_single_tc_class <- function(x) {
  grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R[^]]+\\][ACGT]$", x)
}

# Render one long_visual string as openxlsx2 rich text.
# When brackets=FALSE, bracket characters (<>{}[]) are omitted from the output
# but their contents retain full styling (color, bold, underline).
.rich_long_visual <- function(s, brackets = TRUE, mono_font, cyan) {
  fmt_txt <- openxlsx2::fmt_txt
  if (is.na(s) || !nzchar(s)) {
    return(fmt_txt(""))
  }
  txt <- function(x, ...) fmt_txt(x, font = mono_font, ...)
  br <- function(x) if (brackets) list(txt(x)) else list()
  pieces <- regmatches(
    s,
    gregexpr("<[^>]*>|\\[[^]]*\\]|\\{[^}]*\\}|[^<\\[\\{]+", s, perl = TRUE)
  )[[1]]

  # Pre-pass: derive the repeat unit length from the <...> content by stripping
  # {} delimiters to get the raw sequence length.
  unit_len <- NA_integer_
  for (piece in pieces) {
    if (startsWith(piece, "<")) {
      inner <- substr(piece, 2L, nchar(piece) - 1L)
      unit_len <- nchar(gsub("[{}]", "", inner))
      break
    }
  }

  parts <- list()
  for (p in pieces) {
    if (startsWith(p, "<")) {
      inner <- substr(p, 2L, nchar(p) - 1L)
      parts <- c(parts, list(txt("<")))
      subs <- regmatches(
        inner,
        gregexpr("\\{[^}]*\\}|[^{]+", inner, perl = TRUE)
      )[[1]]
      for (sp in subs) {
        if (startsWith(sp, "{")) {
          parts <- c(
            parts,
            br("{"),
            list(txt(
              substr(sp, 2L, nchar(sp) - 1L),
              bold = TRUE,
              color = cyan,
              underline = "double"
            )),
            br("}")
          )
        } else {
          parts <- c(parts, list(txt(sp, bold = TRUE, color = cyan)))
        }
      }
      parts <- c(parts, list(txt(">")))
    } else if (startsWith(p, "[")) {
      content <- substr(p, 2L, nchar(p) - 1L)
      parts <- c(parts, br("["))
      if (
        !is.na(unit_len) && unit_len > 0L && nchar(content) >= 2L * unit_len
      ) {
        # Split into per-copy chunks; underline even-numbered copies.
        n_copies <- floor(nchar(content) / unit_len)
        for (copy_i in seq_len(n_copies)) {
          copy_seq <- substr(
            content,
            (copy_i - 1L) * unit_len + 1L,
            copy_i * unit_len
          )
          if (copy_i %% 2L == 0L) {
            parts <- c(
              parts,
              list(txt(copy_seq, color = cyan, underline = TRUE))
            )
          } else {
            parts <- c(parts, list(txt(copy_seq, color = cyan)))
          }
        }
        remainder <- substr(content, n_copies * unit_len + 1L, nchar(content))
        if (nzchar(remainder)) {
          parts <- c(parts, list(txt(remainder, color = cyan)))
        }
      } else {
        parts <- c(parts, list(txt(content, color = cyan)))
      }
      parts <- c(parts, br("]"))
    } else if (startsWith(p, "{")) {
      parts <- c(
        parts,
        br("{"),
        list(txt(substr(p, 2L, nchar(p) - 1L), underline = "double")),
        br("}")
      )
    } else {
      parts <- c(parts, list(txt(p)))
    }
  }
  Reduce(`+`, parts)
}

# Write the rosetta table as a formatted Excel workbook via openxlsx2.
.write_rosetta_xlsx <- function(doc, out_path) {
  doc <- data.table::copy(doc)
  wb <- openxlsx2::wb_workbook()$add_worksheet("rosetta")
  wb_dims <- openxlsx2::wb_dims

  # Pre-compute column indices before renaming.
  long_col <- which(names(doc) == "long_visual")
  cosmic_col <- which(names(doc) == "COSMIC_83")
  block_cols <- match(
    c("Koh_476", "Koh_89", "COSMIC_83", "n_indels"),
    names(doc)
  )
  count_col <- which(names(doc) == "n_indels")
  numeric_cols <- which(sapply(doc, is.numeric))

  # Rename display headers.
  data.table::setnames(
    doc,
    old = c(
      "Koh_476",
      "Koh_89",
      "COSMIC_83",
      "n_indels",
      "example_n",
      "long_visual"
    ),
    new = c(
      "476-type",
      "89-type",
      "83-type",
      "N indels",
      "Example n",
      "Indel in context"
    )
  )

  # Write the non-rich columns first (long_col filled in by rich-text pass).
  plain <- data.table::copy(doc)
  plain[[long_col]] <- NA_character_
  wb$add_data(x = plain, na.strings = "")

  # Write long_visual as rich text per cell.
  cyan <- openxlsx2::wb_color(hex = "FFFF1493") # bright pink (DeepPink)
  mono_font <- "Liberation Mono"
  for (i in seq_len(nrow(doc))) {
    wb$add_data(
      x = .rich_long_visual(
        doc[[long_col]][i],
        brackets = FALSE,
        mono_font = mono_font,
        cyan = cyan
      ),
      dims = wb_dims(rows = i + 1L, cols = long_col)
    )
  }

  # Header style: bold, vertically centered; numeric headers rotated vertical.
  wb$add_cell_style(
    dims = wb_dims(rows = 1, cols = seq_along(doc)),
    vertical = "center"
  )
  wb$add_font(dims = wb_dims(rows = 1, cols = seq_along(doc)), bold = "1")
  if (length(numeric_cols) > 0) {
    wb$add_cell_style(
      dims = wb_dims(rows = 1, cols = numeric_cols),
      text_rotation = 90,
      horizontal = "center",
      vertical = "bottom"
    )
  }

  # All body cells: vertically centered; numeric columns horizontally centered.
  body_rows <- seq_len(nrow(doc)) + 1L
  wb$add_cell_style(
    dims = wb_dims(rows = body_rows, cols = seq_along(doc)),
    vertical = "center"
  )
  if (length(numeric_cols) > 0) {
    wb$add_cell_style(
      dims = wb_dims(rows = body_rows, cols = numeric_cols),
      vertical = "center",
      horizontal = "center"
    )
  }

  # Wrap text in 83-type column.
  wb$add_cell_style(
    dims = wb_dims(rows = body_rows, cols = cosmic_col),
    wrap_text = "1",
    vertical = "center"
  )

  # Monospace "Indel in context" cells (font also set per rich-text run, but
  # set the whole-cell font too so plain segments stay aligned).
  wb$add_font(
    dims = wb_dims(rows = body_rows, cols = long_col),
    name = mono_font
  )

  wb$freeze_pane(first_row = TRUE)

  # Vertically merge 476-type / 89-type / 83-type / N indels across each
  # pair's block.
  key <- paste(doc[["476-type"]], doc[["89-type"]], sep = "\r")
  runs <- rle(key)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  for (i in seq_along(runs$lengths)) {
    if (runs$lengths[i] >= 2) {
      rng <- (starts[i]:ends[i]) + 1L
      for (col in block_cols) {
        wb$merge_cells(dims = wb_dims(rows = rng, cols = col))
      }
    }
  }

  # Column widths: numeric columns 4 wide (N indels 7), others by content.
  col_widths <- c(
    17, 18, 16, 7, 3, 70, 12, 12, 10, 8, 8, 12, 8, 8, 8, 10, 8, 8, 8, 8
  )
  col_widths[numeric_cols] <- 4
  col_widths[count_col] <- 7
  wb$set_col_widths(
    cols = seq_along(doc),
    widths = col_widths[seq_along(doc)]
  )

  openxlsx2::wb_save(wb, out_path, overwrite = TRUE)
  invisible(out_path)
}
