# Reverse complement a character vector of DNA strings

Thin vectorized wrapper around
[`Biostrings::reverseComplement()`](https://rdrr.io/pkg/Biostrings/man/reverseComplement.html).
Handles IUPAC ambiguity codes.

## Usage

``` r
revc(x)
```

## Arguments

- x:

  A character vector of DNA sequences.

## Value

A character vector of reverse complements.
