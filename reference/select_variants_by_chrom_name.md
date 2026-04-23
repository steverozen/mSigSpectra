# Restrict a VCF data frame to a user-specified set of chromosome names

Restrict a VCF data frame to a user-specified set of chromosome names

## Usage

``` r
select_variants_by_chrom_name(df, chr.names.to.process, name.of.VCF = NULL)
```

## Arguments

- df:

  An in-memory data frame representing a VCF.

- chr.names.to.process:

  A character vector of chromosome names to keep.

- name.of.VCF:

  Name of the VCF file (for warning messages).

## Value

A list with elements

- `df`: data frame of rows whose `CHROM` is in `chr.names.to.process`.

- `discarded.variants`: **non-NULL only if** any rows were discarded.
