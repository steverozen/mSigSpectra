# Build a one-file Excel table of trimmed indel-context strings, with counts,
# for a chosen set of 476-type mutation classes.
#
# Like build_one_file_rosetta(), this reads one annotated indel VCF (as written
# by annotate_id_vcf()). Unlike it, it looks at *all* rows of the requested
# Koh_476 types (no example sampling), trims each long_visual to left_num
# characters left of "<" and right_num 3' flank characters after the indel
# token ("[...]" repeat text after ">" is dropped, "{...}" microhomology
# text is kept), and
# counts
# each distinct trimmed string.
#
# Usage: devtools::load_all(); source("development/build_476_context_counts.R")

#' @param vcf_path Path to one annotated indel VCF (may be gzipped). Must
#'   contain the columns `Koh_476`, `Koh_89`, and `long_visual`.
#' @param types_476 Character vector of 476-type names to keep, e.g.
#'   `c("Del2:M1", "Del2:U2:R2")`. Closed-interval spellings such as
#'   `"R(5,9)"` are accepted, they are converted to open intervals before
#'   matching.
#' @param out_path Path of the `.xlsx` file to write. If `NULL`, no file is
#'   written and only the table is returned.
#' @param left_num Number of characters to keep to the left of `<`. If `NULL`,
#'   keep the whole 5' flank.
#' @param right_num Number of 3' flank characters to keep after the indel
#'   token. Square-bracketed repeat text after `>` (for example `[AG]`) is
#'   dropped. Curly-braced microhomology text (for example `{A}`) is kept
#'   and does not count toward `right_num`. If `NULL`, keep the whole 3'
#'   flank.
#' @param cap_9 If `TRUE`, drop single-base C/T rows whose repeat count is 10
#'   or more and rewrite `R9` as `R(9,)`, as in build_one_file_rosetta().
#'
#' @return Invisibly, a `data.table` with columns `Koh_476`, `Koh_89`,
#'   `n_indels` (rows in the VCF with that 476-type), `context` (the trimmed
#'   string), and `n` (rows with that trimmed string), sorted in the order of
#'   `types_476`, then descending `n`.
build_476_context_counts <- function(
  vcf_path,
  types_476,
  out_path = NULL,
  left_num = NULL,
  right_num = NULL,
  cap_9 = TRUE
) {
  if (!is.null(out_path) && !requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("Package 'openxlsx2' is required to write the Excel file")
  }
  keep_cols <- c("Koh_476", "Koh_89", "long_visual")

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
    stop(vcf_path, " is missing columns: ", paste(missing, collapse = ", "))
  }
  df$long_visual <- .revcomp_single_tc_long_visual(df$long_visual, df$Koh_476)
  full <- data.table::as.data.table(df[, keep_cols, drop = FALSE])

  if (cap_9) {
    drop <- .is_uncapped_del_repeat(full$Koh_476)
    message("Dropping ", sum(drop), " of ", nrow(full), " uncapped rows")
    full <- full[!drop]
    full$Koh_476 <- sub(
      "((?:Ins|Del)\\([CT]\\)):R9\\]",
      "\\1:R(9,)]",
      full$Koh_476,
      perl = TRUE
    )
  }
  full[, Koh_476 := change_476_type_ids_to_open_intervals(Koh_476)]
  full[, Koh_89 := change_89_type_ids_to_open_intervals(Koh_89)]
  types_476 <- change_476_type_ids_to_open_intervals(types_476)

  not_found <- setdiff(types_476, unique(full$Koh_476))
  if (length(not_found) > 0) {
    warning(
      "476-types not present in ",
      basename(vcf_path),
      ": ",
      paste(not_found, collapse = ", ")
    )
  }
  full <- full[Koh_476 %in% types_476]
  message("Kept ", nrow(full), " rows for ", length(types_476), " 476-types")

  full[, context := .trim_context(long_visual, left_num, right_num)]
  full[, n_indels := .N, by = Koh_476]
  doc <- full[, .(n = .N), by = .(Koh_476, Koh_89, n_indels, context)]

  # Order: the order of types_476, then most common context first.
  doc[, .type_ord := match(Koh_476, types_476)]
  data.table::setorder(doc, .type_ord, Koh_89, -n, context)
  doc[, .type_ord := NULL]

  if (!is.null(out_path)) {
    .write_context_counts_xlsx(doc, out_path)
    message("Wrote ", out_path)
  }
  invisible(doc)
}

# long_visual is "<5'-flank> <indel-token> <3'-flank>", where the indel token
# is "<...>" possibly followed by bracketed repeat or microhomology text such
# as "[AG]" or "{A}". Keep the last left_num characters of the 5' flank, the
# "<...>" token with any "[...]" repeat text after ">" dropped (but "{...}"
# microhomology text kept), and the first
# right_num characters of the 3' flank. Rows without the three-part form are returned with spaces
# removed but otherwise untouched.
.trim_context <- function(x, left_num = NULL, right_num = NULL) {
  out <- gsub(" ", "", x, fixed = TRUE)
  parts <- strsplit(x, " ", fixed = TRUE)
  ok <- !is.na(x) & vapply(parts, length, integer(1)) == 3L
  good <- parts[ok]
  pre <- vapply(good, `[`, character(1), 1)
  mid <- vapply(good, `[`, character(1), 2)
  post <- vapply(good, `[`, character(1), 3)
  # Drop the square-bracketed repeat text after ">", e.g. "[AG]", but keep
  # curly-braced microhomology text such as "{A}".
  mid <- gsub("\\[[^]]*\\]", "", mid)
  if (!is.null(left_num)) {
    pre <- substr(pre, pmax(1, nchar(pre) - left_num + 1), nchar(pre))
  }
  if (!is.null(right_num)) {
    post <- substr(post, 1, right_num)
  }
  out[ok] <- paste0(pre, mid, post)
  out
}

.write_context_counts_xlsx <- function(doc, out_path) {
  doc <- data.table::copy(doc)
  wb <- openxlsx2::wb_workbook()$add_worksheet("contexts")
  wb_dims <- openxlsx2::wb_dims

  ctx_col <- which(names(doc) == "context")
  block_cols <- match(c("Koh_476", "Koh_89", "n_indels"), names(doc))
  numeric_cols <- which(sapply(doc, is.numeric))

  data.table::setnames(
    doc,
    old = c("Koh_476", "Koh_89", "n_indels", "context", "n"),
    new = c("476-type", "89-type", "N indels", "Indel in context", "N")
  )

  plain <- data.table::copy(doc)
  plain[[ctx_col]] <- NA_character_
  wb$add_data(x = plain, na.strings = "")

  cyan <- openxlsx2::wb_color(hex = "FFFF1493")
  mono_font <- "Liberation Mono"
  for (i in seq_len(nrow(doc))) {
    wb$add_data(
      x = .rich_long_visual(
        doc[[ctx_col]][i],
        brackets = TRUE,
        mono_font = mono_font,
        cyan = cyan
      ),
      dims = wb_dims(rows = i + 1L, cols = ctx_col)
    )
  }

  wb$add_font(dims = wb_dims(rows = 1, cols = seq_along(doc)), bold = "1")
  wb$add_cell_style(
    dims = wb_dims(rows = 1, cols = numeric_cols),
    text_rotation = 90,
    horizontal = "center",
    vertical = "bottom"
  )
  body_rows <- seq_len(nrow(doc)) + 1L
  wb$add_cell_style(
    dims = wb_dims(rows = body_rows, cols = seq_along(doc)),
    vertical = "center"
  )
  wb$add_cell_style(
    dims = wb_dims(rows = body_rows, cols = numeric_cols),
    vertical = "center",
    horizontal = "center"
  )
  wb$add_font(dims = wb_dims(rows = body_rows, cols = ctx_col), name = mono_font)
  wb$freeze_pane(first_row = TRUE)

  # Vertically merge 476-type / 89-type / N indels across each block.
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

  col_widths <- c(17, 18, 7, 70, 7)
  wb$set_col_widths(cols = seq_along(doc), widths = col_widths)
  openxlsx2::wb_save(wb, out_path, overwrite = TRUE)
  invisible(out_path)
}
