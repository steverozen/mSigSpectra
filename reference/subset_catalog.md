# Subset a catalog while preserving attributes

Base `[` on a matrix drops attributes other than `dim` / `dimnames`. Use
`subset_catalog()` when you want to keep the catalog's type / ref_genome
/ region / abundance metadata through a subset operation.

## Usage

``` r
subset_catalog(x, rows = NULL, cols = NULL)
```

## Arguments

- x:

  A catalog.

- rows, cols:

  Numeric, logical, or character indices (see
  [base::Extract](https://rdrr.io/r/base/Extract.html)). If `NULL`, all
  rows / columns are kept.

## Value

A catalog (numeric matrix) containing the selected rows and columns,
with the `type`, `counts_or_density`, `ref_genome`, and `region`
attributes of `x` preserved. The `abundance` attribute is preserved only
when all rows are kept, because it is not meaningful for a subset of
mutation types.
