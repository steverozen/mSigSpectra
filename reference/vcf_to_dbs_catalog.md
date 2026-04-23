# Build a DBS mutational-spectrum catalog from an annotated DBS VCF

Build a DBS mutational-spectrum catalog from an annotated DBS VCF

## Usage

``` r
vcf_to_dbs_catalog(
  annotated_vcf,
  type = c("DBS78", "DBS136", "DBS144"),
  ref_genome = NULL,
  region = "unknown",
  sample_name = "count"
)
```

## Arguments

- annotated_vcf:

  A DBS VCF annotated by
  [`annotate_sbs_or_dbs_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_sbs_or_dbs_vcf.md).
  May be the bare annotated `data.table` or the full
  `list(annotated.vcf, discarded.variants)` returned by the annotator.
  Must contain a `seq.<N>bases` column; for `type = "DBS144"` also
  requires `trans.strand` / `bothstrand`.

- type:

  One of `"DBS78"`, `"DBS136"`, `"DBS144"`.

- ref_genome:

  Optional BSgenome object or alias; recorded on the output catalog.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.

- sample_name:

  Column name for the single-sample catalog matrix.

## Value

A single-column numeric matrix with catalog attributes (see
[`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)).
