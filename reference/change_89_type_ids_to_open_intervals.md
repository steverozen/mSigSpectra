# Change 89-type indel category identifiers to use right-open repeat intervals

Replaces bounded repeat-count suffixes like `R(5,9)` with right-open
equivalents like `R(5,)`, making the classification agnostic as to
whether repeats longer than 9 were discarded from the input data.

## Usage

``` r
change_89_type_ids_to_open_intervals(type_89_indel_type_identifiers)
```

## Arguments

- type_89_indel_type_identifiers:

  Character vector of 89-type indel category identifiers.

## Value

Character vector the same length as the input, with bounded repeat
suffixes replaced by right-open equivalents.

## Examples

``` r
change_89_type_ids_to_open_intervals(
  c("Ins(2,):R(5,9)", "Ins(C):R(7,9)", "[Del(T):R(8,9)]",
    "[Ins(T):R(8,9)]", "Del(2,):U(1,2):R(5,9)", "Del(3,):U(3,):R(3,9)")
)
#> [1] "Ins(2,):R(5,)"        "Ins(C):R(7,)"         "[Del(T):R(8,)]"      
#> [4] "[Ins(T):R(8,)]"       "Del(2,):U(1,2):R(5,)" "Del(3,):U(3,):R(3,)" 
```
