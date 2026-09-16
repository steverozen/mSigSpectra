# Check whether an object looks like an mSigSpectra catalog

Returns `TRUE` if `x` is a numeric matrix with the five catalog
attributes (`type`, `counts_or_density`, `ref_genome`, `region`,
`abundance`) and canonical rownames for its type.

## Usage

``` r
is_catalog(x)
```

## Arguments

- x:

  Any R object.

## Value

A single logical value: `TRUE` if `x` is a numeric matrix carrying the
catalog attributes with the canonical row names for its `type`,
otherwise `FALSE`.
