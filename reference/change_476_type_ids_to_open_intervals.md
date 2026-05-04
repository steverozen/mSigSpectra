# Change 476-type indel category identifiers to use right-open repeat intervals

Replaces bounded repeat-count suffixes like `:R(5,9)` with right-open
equivalents like `:R(5,)`, making the classification agnostic as to
whether repeats longer than 9 were discarded from the input data.

## Usage

``` r
change_476_type_ids_to_open_intervals(type_476_indel_type_identifiers)
```

## Arguments

- type_476_indel_type_identifiers:

  Character vector of 476-type indel category identifiers, e.g. as
  returned by `categorize_indels_in_vcf()`.

## Value

Character vector the same length as the input, with `:R(5,9)` at the end
of each string replaced by `:R(5,)`.

## Examples

``` r
change_476_type_ids_to_open_intervals(
  c("Del(C):Ins(C):R(5,9)", "Del(T):R(3,5)", "Ins(C):R(5,9)")
)
#> [1] "Del(C):Ins(C):R(5,)" "Del(T):R(3,5)"       "Ins(C):R(5,)"       
```
