# Read transcript ranges from a GENCODE-derived CSV

Reads a CSV with columns `chrom`, `start`, `end`, `strand`,
`Ensembl.gene.ID`, `gene.symbol` (1-based coordinates), orders the
chromosome factor canonically (1..22/19, X, Y), and returns a keyed
`data.table`.

## Usage

``` r
read_transcript_ranges(file)
```

## Arguments

- file:

  Path to the transcript-range CSV file.

## Value

A `data.table` keyed on `chrom`, `start`, `end`.
