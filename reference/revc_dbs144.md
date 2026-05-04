# Reverse complement strings that represent stranded DBSs

Input is a 4-character string where characters 1-2 are the REF
dinucleotide and characters 3-4 are the ALT (e.g. `"AATC"` = `AA>TC`).
Returns the reverse complement of the first 2 characters concatenated
with the reverse complement of the last 2 characters, e.g. `"AATC"`
returns `"TTGA"`.

## Usage

``` r
revc_dbs144(mutstring)
```

## Arguments

- mutstring:

  A character vector of 4-letter strings.
