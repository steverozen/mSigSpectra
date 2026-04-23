# Combine catalogs across samples (column-bind)

Checks that all input catalogs share the same `type`,
`counts_or_density`, `ref_genome`, and `region`, then `cbind`s their
matrices and re-applies the shared attributes.

## Usage

``` r
cbind_catalogs(catalogs)
```

## Arguments

- catalogs:

  A list of catalogs.

## Value

A single catalog with `ncol` equal to the sum of the input `ncol`s.
