# Compute length of longest common prefix of two strings

Fast alternative to
[`Biostrings::lcprefix`](https://rdrr.io/pkg/Biostrings/man/lcsuffix.html)
that avoids S4 method dispatch overhead.

## Usage

``` r
lcprefix_fast(a, b)
```

## Arguments

- a:

  A single character string.

- b:

  A single character string.

## Value

An integer: the number of leading characters that match.
