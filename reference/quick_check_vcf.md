# Filter, deduplicate, and check an annotated indel VCF

Helper used by
[`annot_vcf_to_83_catalog`](https://steverozen.github.io/mSigSpectra/reference/annot_vcf_to_83_catalog.md),
[`annot_vcf_to_89_catalog`](https://steverozen.github.io/mSigSpectra/reference/annot_vcf_to_89_catalog.md),
and
[`annot_vcf_to_476_catalog`](https://steverozen.github.io/mSigSpectra/reference/annot_vcf_to_476_catalog.md)
to perform shared preprocessing: rename `#CHROM` to `CHROM`, optionally
filter to PASS variants, warn about positions with differing ALT
alleles, and deduplicate by position.

## Usage

``` r
quick_check_vcf(annot_vcf, FILTER_PASS = FALSE, do_message = FALSE)
```

## Arguments

- annot_vcf:

  A data frame with at least columns `CHROM` (or `#CHROM`), `POS`, and
  `ALT`. If `FILTER_PASS` is `TRUE`, a `FILTER` column is also required.

- FILTER_PASS:

  If `TRUE`, retain only rows where the `FILTER` column equals `"PASS"`.

- do_message:

  If `TRUE`, emit diagnostic messages showing row counts at each
  processing step.

## Value

A data frame deduplicated by position, with a `pos_id` column added.
