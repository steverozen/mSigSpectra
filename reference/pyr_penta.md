# Normalize SBS1536 pentanucleotide + ALT strings to pyrimidine form

Input strings are 6 characters: a 5-base pentanucleotide context
followed by a 1-base ALT (e.g. `"ATGCTT"` = `ATGCT>T`). If the center of
the pentanucleotide (position 3) is `A` or `G`, the pentanucleotide and
the ALT are reverse-complemented so that the center becomes `C` or `T` —
the pyrimidine form used in canonical SBS96/1536 row names.

## Usage

``` r
pyr_penta(x)
```

## Arguments

- x:

  A character vector of 6-letter strings.
