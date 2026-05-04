# Reverse complement strings that represent stranded SBSs

Input is a 4-character string where characters 1-3 are the trinucleotide
context and character 4 is the ALT (e.g. `"AATC"` = `AAT>ACT`). Returns
the reverse complement of the first 3 characters concatenated with the
reverse complement of the last character, e.g. `"AATC"` returns
`"ATTG"`.

## Usage

``` r
revc_sbs96(mutstring)
```

## Arguments

- mutstring:

  A character vector of 4-letter strings.
