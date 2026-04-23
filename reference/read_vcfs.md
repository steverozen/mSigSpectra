# Read multiple VCF files

Thin batch wrapper around
[`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md).

## Usage

``` r
read_vcfs(files, ..., names_of_vcfs = NULL)
```

## Arguments

- files:

  Character vector of file paths / URLs.

- ...:

  Passed through to
  [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md).

- names_of_vcfs:

  Optional character vector of names (same length as `files`). Defaults
  to stripping extensions from basenames.

## Value

A named list of `data.table`s, one per input file.
