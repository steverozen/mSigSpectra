# Collapse a higher-resolution catalog to a lower-resolution one

Supports:

- `SBS1536 -> SBS96`: collapse pentanucleotide → trinucleotide contexts.

- `SBS192 -> SBS96`: sum across the two transcript strands.

- `DBS144 -> DBS78`: sum across the two transcript strands.

## Usage

``` r
collapse_catalog(catalog, to = c("SBS96", "DBS78"))
```

## Arguments

- catalog:

  An mSigSpectra catalog.

- to:

  Target catalog type (`"SBS96"`, `"DBS78"`).

## Value

A new catalog with the collapsed rows.
