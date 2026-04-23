# Validate a `region` argument

Validate a `region` argument

## Usage

``` r
stop_if_region_illegal(region)
```

## Arguments

- region:

  Character string to check.

## Value

`NULL` invisibly; raises an error if `region` is not one of `"genome"`,
`"exome"`, `"transcript"`, `"unknown"`.
