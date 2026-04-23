# Write a mutational-spectrum catalog to a file

Writes `catalog` in ICAMS-native CSV format (the only format currently
supported for writing). The row-header columns that precede the sample
columns are taken from the shipped `catalog.row.headers` object for the
catalog's type.

## Usage

``` r
write_catalog(
  catalog,
  file,
  format = c("ICAMS", "SigProfiler", "COSMIC"),
  sep = ","
)
```

## Arguments

- catalog:

  A catalog (see
  [`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)).

- file:

  Output path.

- format:

  Output format. Currently only `"ICAMS"` is supported; `"SigProfiler"`
  and `"COSMIC"` error with an informative message.

- sep:

  Column separator. Defaults to `","` for ICAMS format.

## Value

`file` invisibly.
