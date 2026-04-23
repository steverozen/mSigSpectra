# Build an SBS mutational-spectrum catalog from an annotated SBS VCF

Returns a single catalog matrix of the requested `type`. Intermediate
matrices for the other SBS resolutions are still computed (cheap) but
not returned, keeping the public API focused on "one call, one catalog
type".

## Usage

``` r
vcf_to_sbs_catalog(
  annotated_vcf,
  type = c("SBS96", "SBS192", "SBS1536"),
  ref_genome = NULL,
  region = "unknown",
  sample_name = "count"
)
```

## Arguments

- annotated_vcf:

  An SBS VCF annotated by
  [`annotate_sbs_or_dbs_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_sbs_or_dbs_vcf.md).
  May be the bare annotated `data.table` or the full
  `list(annotated.vcf, discarded.variants)` returned by the annotator.
  Must contain a `seq.<N>bases` column; for `type = "SBS192"` also
  requires `trans.strand` / `bothstrand`.

- type:

  One of `"SBS96"`, `"SBS192"`, `"SBS1536"`.

- ref_genome:

  Optional BSgenome object or alias; recorded on the output catalog.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.

- sample_name:

  Column name for the single-sample catalog matrix.

## Value

A single-column numeric matrix with catalog attributes (see
[`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)).
