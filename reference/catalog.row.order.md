# Canonical row orders for catalogs

Row-name orders for every catalog type supported by mSigSpectra. Exposed
for users who need to build catalogs from formats not otherwise
supported.

## Usage

``` r
catalog.row.order
```

## Format

A named list of character vectors keyed by catalog type (e.g. `"SBS96"`,
`"DBS78"`, `"ID83"`).

## Examples

``` r
catalog.row.order$SBS96[1:5]
#> [1] "ACAA" "ACCA" "ACGA" "ACTA" "CCAA"
```
