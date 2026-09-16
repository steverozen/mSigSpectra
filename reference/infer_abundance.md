# Infer k-mer abundance from attributes

Looks up the appropriate entry of `all.abundance` keyed on `ref_genome`
and `region`. Returns `NULL` if no entry is found, if the catalog is a
signature / density catalog with a flat abundance, or if the catalog is
COMPOSITE (which has no meaningful abundance).

## Usage

``` r
infer_abundance(x, ref_genome, region, counts_or_density)
```

## Value

A named integer vector of k-mer counts (the abundance for the catalog's
context size), or `NULL` if no abundance applies.
