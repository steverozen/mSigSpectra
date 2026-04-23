# Standardize chromosome names in the first column of a data frame

Drops rows whose chromosome name contains any of `GL`, `KI`, `random`,
`Hs` (anchored at start), `M`, or `JH`, and strips any leading `chr`
prefix from the remainder.

## Usage

``` r
standard_chrom_name(df)
```

## Arguments

- df:

  A data frame whose first column contains chromosome names.

## Value

`df` restricted to rows with canonical chromosome names (`1:22`, `X`,
`Y`), with `chr` prefixes removed.
