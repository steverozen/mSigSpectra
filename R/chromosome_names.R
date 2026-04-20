#' Standardize chromosome names in the first column of a data frame
#'
#' Drops rows whose chromosome name contains any of `GL`, `KI`, `random`,
#' `Hs` (anchored at start), `M`, or `JH`, and strips any leading `chr`
#' prefix from the remainder.
#'
#' @param df A data frame whose first column contains chromosome names.
#'
#' @return `df` restricted to rows with canonical chromosome names
#'   (`1:22`, `X`, `Y`), with `chr` prefixes removed.
#'
#' @keywords internal
standard_chrom_name <- function(df) {
  patterns <- c("GL", "KI", "random", "^Hs", "M", "JH")
  for (p in patterns) {
    hits <- grepl(p, df[[1]])
    if (any(hits)) df <- df[!hits, ]
  }
  df[[1]] <- sub(pattern = "chr", replacement = "", df[[1]])
  df
}

#' Standardize the chromosome names in a VCF data.frame
#'
#' Splits `df` into rows with canonical chromosome names and rows with
#' non-standard chromosome names (those containing any of `GL`, `KI`,
#' `random`, `Hs`, `M`, `JH`, `fix`, `alt`). Emits a warning when any
#' rows are discarded.
#'
#' @param df An in-memory data frame representing a VCF; must contain a
#'   `CHROM` column.
#' @param name.of.VCF Name of the VCF file (for warning messages).
#'
#' @return A list with elements
#'   * `df`: data frame of rows with canonical chromosome names.
#'   * `discarded.variants`: **non-NULL only if** any rows were discarded;
#'     each discarded row gets a `discarded.reason` column.
#' @md
#'
#' @keywords internal
standard_chrom_name_new <- function(df, name.of.VCF = NULL) {
  discarded.variants <- df[0, ]
  patterns <- list(
    list(pat = "GL",   label = "GL"),
    list(pat = "KI",   label = "KI"),
    list(pat = "random", label = "random"),
    list(pat = "^Hs",  label = "Hs"),
    list(pat = "M",    label = "M"),
    list(pat = "JH",   label = "JH"),
    list(pat = "fix",  label = "fix"),
    list(pat = "alt",  label = "alt")
  )

  n.orig <- nrow(df)
  for (p in patterns) {
    hits <- grepl(p$pat, df$CHROM)
    if (any(hits)) {
      warning(
        "In VCF ",
        ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
        " ", sum(hits), " row out of ", n.orig,
        " had chromosome names that contain '", p$label, "' and were removed. ",
        "See discarded.variants in the return value for more details"
      )
      to.remove <- df[hits, ]
      to.remove$discarded.reason <-
        paste0("Chromosome name contains \"", p$label, "\"")
      discarded.variants <- rbind(discarded.variants, to.remove)
      df <- df[!hits, ]
    }
  }

  if (nrow(discarded.variants) == 0) {
    list(df = df)
  } else {
    list(df = df, discarded.variants = discarded.variants)
  }
}

#' Restrict a VCF data frame to a user-specified set of chromosome names
#'
#' @param df An in-memory data frame representing a VCF.
#' @param chr.names.to.process A character vector of chromosome names to keep.
#' @param name.of.VCF Name of the VCF file (for warning messages).
#'
#' @return A list with elements
#'   * `df`: data frame of rows whose `CHROM` is in `chr.names.to.process`.
#'   * `discarded.variants`: **non-NULL only if** any rows were discarded.
#' @md
#'
#' @keywords internal
select_variants_by_chrom_name <-
  function(df, chr.names.to.process, name.of.VCF = NULL) {
    keep <- df$CHROM %in% chr.names.to.process
    df1 <- df[keep, ]
    discarded.variants <- df[!keep, ]

    if (nrow(discarded.variants) == 0) {
      list(df = df1)
    } else {
      warning(
        "In VCF ",
        ifelse(is.null(name.of.VCF), "", dQuote(name.of.VCF)),
        " ", nrow(discarded.variants), " row out of ", nrow(df),
        " had chromosome names that were not selected by user and were removed. ",
        "See discarded.variants in the return value for more details"
      )
      discarded.variants$discarded.reason <-
        "Chromosome names not selected by user"
      list(df = df1, discarded.variants = discarded.variants)
    }
  }

#' Check and, if possible, correct the chromosome names in a VCF data.frame
#'
#' Harmonizes the `CHROM` column of `vcf.df` with the chromosome names used
#' by `ref.genome` (a BSgenome object). Adds or strips the `chr` prefix as
#' needed; for human and mouse, translates `23/24` or `20/21` to `X/Y`.
#'
#' @param vcf.df A VCF as a data frame with a `CHROM` column.
#' @param ref.genome A BSgenome object (e.g. `BSgenome.Hsapiens.UCSC.hg38`).
#' @param name.of.VCF Name of the VCF file (for warning/error messages).
#'
#' @return A character vector of chromosome names that can be used as a
#'   replacement for `vcf.df$CHROM`. Errors if reconciliation is impossible.
#'
#' @keywords internal
check_and_fix_chrom_names <- function(vcf.df, ref.genome, name.of.VCF = NULL) {
  names.to.check <- unique(vcf.df$CHROM)

  if (!sum(grepl("^chr", names.to.check)) %in% c(0, length(names.to.check))) {
    stop(
      "\nNaming of chromosomes in VCF ", dQuote(name.of.VCF),
      " is not consistent: ", paste(names.to.check, collapse = " ")
    )
  }

  ref.genome.names <- GenomeInfoDb::seqnames(ref.genome)

  not.matched <- setdiff(names.to.check, ref.genome.names)
  if (length(not.matched) == 0) return(vcf.df$CHROM)

  vcf.has.chr.prefix <- any(grepl("^chr", names.to.check))
  ref.has.chr.prefix <- any(grepl("^chr", ref.genome.names))

  new.chr.names <- vcf.df$CHROM
  if (ref.has.chr.prefix && !vcf.has.chr.prefix) {
    names.to.check <- paste0("chr", names.to.check)
    new.chr.names <- paste0("chr", new.chr.names)
    if (length(setdiff(names.to.check, ref.genome.names)) == 0) {
      return(new.chr.names)
    }
  }

  if (!ref.has.chr.prefix && vcf.has.chr.prefix) {
    names.to.check <- gsub("chr", "", names.to.check)
    new.chr.names <- gsub("chr", "", new.chr.names)
    if (length(setdiff(names.to.check, ref.genome.names)) == 0) {
      return(new.chr.names)
    }
  }

  organism <- BSgenome::organism(ref.genome)

  check_for_possible_match <- function(chr1, chr2) {
    if (chr1 %in% names.to.check) {
      if (chr2 %in% names.to.check) {
        emit_warning <- function(x, y) {
          warning(
            "\n", x, " and ", y, " both are chromosome names in VCF ",
            dQuote(name.of.VCF),
            "for ", organism, ". ", x, " has been changed to ", y,
            " internally for downstream processing"
          )
        }
        if (vcf.has.chr.prefix) {
          if (grepl("^chr", chr1)) {
            emit_warning(chr1, chr2)
          } else {
            emit_warning(paste0("chr", chr1), paste0("chr", chr2))
          }
        } else {
          if (!grepl("^chr", chr1)) {
            emit_warning(chr1, chr2)
          } else {
            emit_warning(gsub("chr", "", chr1), gsub("chr", "", chr2))
          }
        }
      }
      new.chr.names[new.chr.names == chr1] <<- chr2
      names.to.check <- setdiff(names.to.check, chr1)
      names.to.check <<- unique(c(names.to.check, chr2))
    }
  }

  if (organism == "Homo sapiens") {
    check_for_possible_match("chr23", "chrX")
    check_for_possible_match("chr24", "chrY")
    check_for_possible_match("23", "X")
    check_for_possible_match("24", "Y")
  }
  if (organism == "Mus musculus") {
    check_for_possible_match("chr20", "chrX")
    check_for_possible_match("chr21", "chrY")
    check_for_possible_match("20", "X")
    check_for_possible_match("21", "Y")
  }

  if (length(setdiff(names.to.check, ref.genome.names)) == 0) {
    return(new.chr.names)
  }

  stop(
    "\nChromosome names in VCF ", dQuote(name.of.VCF),
    " not in ref.genome for ", organism, ": ",
    paste(not.matched, collapse = " ")
  )
}

#' Check and, if possible, correct the chromosome names in a trans.ranges table
#'
#' Harmonizes `trans.ranges$chrom` with the chromosome naming used by
#' `vcf.df$CHROM`, adding or stripping `chr` as needed and translating
#' organism-specific numeric X/Y encodings.
#'
#' @param trans.ranges A `data.table` of transcript ranges (see
#'   [trans.ranges]).
#' @param vcf.df A VCF as a data frame with a `CHROM` column.
#' @param ref.genome A BSgenome object used to determine organism.
#' @param name.of.VCF Name of the VCF file.
#'
#' @return A character vector of chromosome names that can be used as a
#'   replacement for `trans.ranges$chrom`. Errors if reconciliation is
#'   impossible.
#'
#' @keywords internal
check_and_fix_chrom_names_for_trans_ranges <-
  function(trans.ranges, vcf.df, ref.genome, name.of.VCF = NULL) {
    names.to.check <- as.character(unique(trans.ranges$chrom))

    vcf.chr.names <- as.character(unique(vcf.df$CHROM))
    if (!sum(grepl("^chr", vcf.chr.names)) %in% c(0, length(vcf.chr.names))) {
      stop(
        "\nNaming of chromosomes in VCF ", dQuote(name.of.VCF),
        " is not consistent: ", paste(vcf.chr.names, collapse = " ")
      )
    }

    not.matched <- setdiff(vcf.chr.names, names.to.check)
    if (length(not.matched) == 0) return(trans.ranges$chrom)

    vcf.has.chr.prefix <- any(grepl("^chr", vcf.chr.names))
    trans.has.chr.prefix <- any(grepl("^chr", names.to.check))

    new.chr.names <- as.character(trans.ranges$chrom)
    if (trans.has.chr.prefix && !vcf.has.chr.prefix) {
      names.to.check <- gsub("chr", "", names.to.check)
      new.chr.names <- gsub("chr", "", names.to.check)
      names.to.check <- paste0("chr", names.to.check)
      new.chr.names <- paste0("chr", new.chr.names)
      if (length(setdiff(vcf.chr.names, names.to.check)) == 0) {
        return(new.chr.names)
      }
    }

    if (!trans.has.chr.prefix && vcf.has.chr.prefix) {
      names.to.check <- paste0("chr", names.to.check)
      new.chr.names <- paste0("chr", new.chr.names)
      if (length(setdiff(vcf.chr.names, names.to.check)) == 0) {
        return(new.chr.names)
      }
    }

    organism <- BSgenome::organism(ref.genome)

    check_for_possible_match <- function(chr1, chr2) {
      if (chr2 %in% vcf.chr.names) {
        if (chr1 %in% vcf.chr.names) {
          emit_stop <- function() {
            stop(
              "\n", chr2, " and ", chr1, " both are chromosome names in VCF ",
              dQuote(name.of.VCF),
              ", which should not be the case for ", organism,
              ". Please check your data or specify the correct ref.genome argument"
            )
          }
          if (vcf.has.chr.prefix) {
            if (grepl("^chr", chr2)) emit_stop() else emit_stop()
          } else {
            if (!grepl("^chr", chr1)) emit_stop() else emit_stop()
          }
        }
        new.chr.names[new.chr.names == chr1] <<- chr2
        names.to.check <- setdiff(names.to.check, chr1)
        names.to.check <<- unique(c(names.to.check, chr2))
      }
    }

    if (organism == "Homo sapiens") {
      check_for_possible_match("chrX", "chr23")
      check_for_possible_match("chrY", "chr24")
      check_for_possible_match("X", "23")
      check_for_possible_match("Y", "24")
    }
    if (organism == "Mus musculus") {
      check_for_possible_match("chrX", "chr20")
      check_for_possible_match("chrY", "chr21")
      check_for_possible_match("X", "20")
      check_for_possible_match("Y", "21")
    }

    if (length(setdiff(vcf.chr.names, names.to.check)) == 0) {
      return(new.chr.names)
    }

    stop(
      "\nChromosome names in VCF ", dQuote(name.of.VCF),
      " not in trans ranges for ", organism, ": ",
      paste(not.matched, collapse = " ")
    )
  }
