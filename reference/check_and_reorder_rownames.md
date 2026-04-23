# Reorder catalog rows to the canonical order for its type

Reorder catalog rows to the canonical order for its type

## Usage

``` r
check_and_reorder_rownames(x, type)
```

## Arguments

- x:

  A matrix with rownames.

- type:

  Catalog type (e.g. `"SBS96"`).

## Value

`x` with rows reordered to match `catalog.row.order[[type]]`. Errors if
any canonical rowname is missing.
