# Build an ID (indel) mutational-spectrum catalog from an annotated ID VCF

Turns an indel-annotated VCF (with `COSMIC_83` / `Koh_89` / `Koh_476`
columns as produced by
[`annotate_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_id_vcf.md))
into a count matrix for the requested ID classification scheme.

## Usage

``` r
vcf_to_id_catalog(
  annotated_vcf,
  type = c("ID83", "ID89", "ID476"),
  ref_genome = NULL,
  region = "unknown",
  sample_name = "count",
  FILTER_PASS = TRUE,
  clip_le_9 = TRUE
)
```

## Arguments

- annotated_vcf:

  An ID VCF annotated by
  [`annotate_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_id_vcf.md).
  May be the bare annotated `data.table` or the full
  `list(annotated.vcf, discarded.variants)` returned by the annotator.
  Must contain the categorization column corresponding to `type`.

- type:

  One of `"ID83"`, `"ID89"`, `"ID476"`.

- ref_genome:

  Optional BSgenome object or alias; recorded on the output catalog.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.

- sample_name:

  Column name for the single-sample catalog matrix.

- FILTER_PASS:

  If `TRUE`, retain only rows where the VCF `FILTER` column is `"PASS"`.

- clip_le_9:

  If `TRUE`, drop variants with repeat count `R > 9`, approximating
  PCAWG indel calling.

## Value

A single-column numeric matrix with catalog attributes (see
[`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)).
