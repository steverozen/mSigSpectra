# Infer catalog type from the number of rows of a matrix

Infer catalog type from the number of rows of a matrix

## Usage

``` r
infer_catalog_type(n_rows)
```

## Arguments

- n_rows:

  Number of rows.

## Value

A type identifier (e.g. `"SBS96"`, `"DBS78"`, `"ID83"`). Errors if
`n_rows` does not correspond to a supported catalog type.
