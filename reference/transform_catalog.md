# Transform a catalog between counts and density

Converts between counts (raw mutation counts per category) and density
(mutations per megabase of context) representations, and between catalog
↔ signature (column-normalized) forms. The transformation is expressed
by multiplying each category's count by
`target_abundance[source.n.mer] / source_abundance[source.n.mer]`, where
`source.n.mer` is the reference context encoded in the row name (first 3
characters for SBS96/192, 5 for SBS1536, 2 for DBS78/144, 4 for DBS136).

## Usage

``` r
transform_catalog(
  catalog,
  target_ref_genome = NULL,
  target_region = NULL,
  target_counts_or_density = NULL,
  target_abundance = NULL
)
```

## Arguments

- catalog:

  An mSigSpectra catalog (see
  [`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md))
  with attributes `type`, `counts_or_density`, `ref_genome`, `region`,
  `abundance`.

- target_ref_genome, target_region, target_counts_or_density:

  Target attributes. If `NULL`, the catalog's current attribute is used.

- target_abundance:

  Optional target abundance vector. If `NULL`, inferred from
  `target_ref_genome` + `target_region`.

## Value

A new catalog (plain matrix with attributes).
