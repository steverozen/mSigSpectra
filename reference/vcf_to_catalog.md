# Build a mutational-spectrum catalog from an annotated VCF

Dispatches to a type-specific internal builder based on `type`. Supports
SBS96 / SBS192 / SBS1536 today; DBS and ID catalog types are stubs
pending the DBS- and indel-specific ports.

## Usage

``` r
vcf_to_catalog(
  annotated_vcf,
  type = c("SBS96", "SBS192", "SBS1536", "DBS78", "DBS136", "DBS144", "ID83", "ID89",
    "ID166", "ID476"),
  ref_genome = NULL,
  region = "unknown",
  sample_name = "count"
)
```

## Arguments

- annotated_vcf:

  A VCF annotated by
  [`annotate_sbs_or_dbs_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_sbs_or_dbs_vcf.md)
  (for SBS / DBS types) or
  [`annotate_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_id_vcf.md)
  (for ID types). May be the bare annotated `data.table` or the full
  `list(annotated.vcf, discarded.variants)` returned by either
  annotator. For SBS types the table must contain a `seq.<N>bases`
  column and (for SBS192) `trans.strand` / `bothstrand`.

- type:

  Catalog type identifier (one of `"SBS96"`, `"SBS192"`, `"SBS1536"`,
  `"DBS78"`, `"DBS136"`, `"DBS144"`, `"ID83"`, `"ID89"`, `"ID166"`,
  `"ID476"`).

- ref_genome:

  Optional BSgenome object or alias; recorded on the output catalog.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.

- sample_name:

  Column name for the single-sample catalog matrix (default `"count"`).

## Value

A single-column numeric matrix with catalog attributes (see
[`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)).
