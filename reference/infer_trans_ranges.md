# Infer transcript ranges for a reference genome

If `trans_ranges` is supplied, it is returned unchanged. Otherwise,
returns the shipped transcript-ranges table for the given reference
genome (GRCh37 / GRCh38 / GRCm38).

## Usage

``` r
infer_trans_ranges(ref_genome, trans_ranges = NULL)
```

## Arguments

- ref_genome:

  A BSgenome object or a character identifier accepted by
  [`normalize_genome_arg()`](https://steverozen.github.io/mSigSpectra/reference/normalize_genome_arg.md).

- trans_ranges:

  Optional user-supplied transcript ranges.

## Value

A `data.table` of transcript ranges, or `NULL` if no shipped table is
available and none was supplied.
