# Standardize the chromosome names in a VCF data.frame

Splits `df` into rows with canonical chromosome names and rows with
non-standard chromosome names (those containing any of `GL`, `KI`,
`random`, `Hs`, `M`, `JH`, `fix`, `alt`). Emits a warning when any rows
are discarded.

## Usage

``` r
standard_chrom_name_new(df, name.of.VCF = NULL)
```

## Arguments

- df:

  An in-memory data frame representing a VCF; must contain a `CHROM`
  column.

- name.of.VCF:

  Name of the VCF file (for warning messages).

## Value

A list with elements

- `df`: data frame of rows with canonical chromosome names.

- `discarded.variants`: **non-NULL only if** any rows were discarded;
  each discarded row gets a `discarded.reason` column.
