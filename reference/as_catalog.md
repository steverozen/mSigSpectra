# Turn a numeric matrix into a mutational-spectrum catalog

Attaches the standard catalog attributes (`type`, `counts_or_density`,
`ref_genome`, `region`, `abundance`) to `x` and returns it. **No S3
class is set** — mSigSpectra catalogs are plain matrices with
attributes; functions key off `attr(x, "type")` via explicit checks.

## Usage

``` r
as_catalog(
  x,
  type = NULL,
  ref_genome = NULL,
  region = "unknown",
  abundance = NULL,
  counts_or_density = "counts",
  infer_rownames = FALSE
)
```

## Arguments

- x:

  A numeric matrix, data.frame coercible to numeric, or a named numeric
  vector (converted to a one-column matrix with the vector names as
  rownames).

- type:

  Optional catalog type identifier (`"SBS96"`, `"DBS78"`, `"ID83"`,
  etc.). Inferred from `nrow(x)` if `NULL`.

- ref_genome:

  Optional BSgenome object or alias; recorded as an attribute and used
  for abundance lookup.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`. For
  stranded catalogs (`SBS192`, `DBS144`) `"genome"` is silently promoted
  to `"transcript"`.

- abundance:

  Optional named numeric vector of k-mer counts. If `NULL`, inferred
  from shipped `all.abundance` keyed on `ref_genome` and `region`.

- counts_or_density:

  One of `"counts"`, `"density"`, `"counts.signature"`,
  `"density.signature"`.

- infer_rownames:

  If `TRUE` and `x` has no rownames, the canonical rownames for the
  catalog type are attached (**assuming the row order is already
  correct**). If `FALSE`, `x` must already have the canonical rownames.

## Value

`x` as a numeric matrix with attributes set.

## Examples

``` r
m <- matrix(
  1, nrow = 96, ncol = 1,
  dimnames = list(catalog.row.order$SBS96, "sample1")
)
cat96 <- as_catalog(m)
attr(cat96, "type")   # "SBS96"
#> [1] "SBS96"
attr(cat96, "region") # "unknown"
#> [1] "unknown"
```
