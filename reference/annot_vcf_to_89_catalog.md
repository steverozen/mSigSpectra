# Convert an annotated indel VCF to a Koh 89-category catalog

Take an annotated indel VCF data frame (with column `Koh_89` as produced
by indel classification functions) and produce a single-column data
frame of mutation counts in the 89-category Koh classification scheme.

## Usage

``` r
annot_vcf_to_89_catalog(
  annot_vcf,
  sample_id = "no_sample_id_provided",
  FILTER_PASS = FALSE,
  do_message = FALSE,
  clip_le_9 = FALSE
)
```

## Arguments

- annot_vcf:

  A data frame with at least columns `CHROM`, `POS`, `ALT`, and
  `Koh_89`. If `FILTER_PASS` is `TRUE`, a `FILTER` column is also
  required.

- sample_id:

  A character string used as the column name in the returned data frame.

- FILTER_PASS:

  If `TRUE`, retain only rows where the `FILTER` column equals `"PASS"`.

- do_message:

  If `TRUE`, emit diagnostic messages showing row counts at each
  processing step.

- clip_le_9:

  Only keep variants with "R" \<= 9, to approximate PCAWG indel calling.

## Value

A single-column data frame with 89 rows (one per Koh category) and
integer mutation counts. Row names are the Koh 89 category strings; the
column name is `sample_id`.

## Details

The function:

1.  Optionally filters to PASS variants.

2.  Removes duplicate positions (warns if ALT alleles differ).

3.  Tallies counts per Koh 89 category and returns a data frame with one
    row per category (using `mSigSpectra::catalog.row.order$ID89`).
