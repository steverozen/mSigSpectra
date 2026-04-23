# Validate a `region` argument for catalog types that require transcript strand

SBS192, DBS144 and similar stranded catalogs cannot be built from
`region = "genome"`, since variants outside transcripts have no strand.

## Usage

``` r
stop_if_transcribed_region_illegal(region)
```

## Arguments

- region:

  Character string to check.
