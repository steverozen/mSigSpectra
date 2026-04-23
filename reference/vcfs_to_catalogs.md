# Convenience wrapper: read, annotate, and catalog a set of VCF files

Thin composition of
[`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md) +
[`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md) +
[`annotate_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_vcf.md) +
[`vcf_to_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_catalog.md)
for a list of files.

## Usage

``` r
vcfs_to_catalogs(
  files,
  types,
  ref_genome,
  region = "genome",
  trans_ranges = NULL,
  names_of_vcfs = NULL
)
```

## Arguments

- files:

  Character vector of VCF paths / URLs.

- types:

  Character vector of catalog types to produce (e.g.
  `c("SBS96", "SBS192")`). Each file contributes one column per type.

- ref_genome:

  A BSgenome object or character alias.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.

- trans_ranges:

  Optional transcript ranges (see
  [`add_transcript_strand()`](https://steverozen.github.io/mSigSpectra/reference/add_transcript_strand.md)).

- names_of_vcfs:

  Optional character vector of names (same length as `files`); defaults
  to basenames without extension.

## Value

A named list keyed by `types`. Each element is a multi-sample catalog
matrix (columns = input files).
