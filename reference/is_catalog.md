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
