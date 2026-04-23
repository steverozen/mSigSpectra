# Segment a single indel sequence using Rcpp interface

Segment a single indel sequence using Rcpp interface

## Usage

``` r
seg_simple(ins_or_del, string, context)
```

## Arguments

- ins_or_del:

  Character indicating insertion ("i") or deletion ("d")

- string:

  A single indel sequence to segment (character scalar)

- context:

  A single flanking context sequence (character scalar)

## Value

A list with the following elements:

- unit:

  Repeat unit sequence (character)

- unit_length:

  Unit length (integer)

- internal_rep:

  Internal repeat region (character)

- internal_reps:

  Internal repeat count (integer)

- spacer:

  Spacer sequence (character)

- spacer_length:

  Spacer length (integer)

- prime3_rep:

  3' flanking repeat region (character)

- prime3_reps:

  3' flanking repeat count (integer)

- original_reps:

  Original repeat count (integer)

## Details

This function segments an indel sequence by finding the optimal repeat
unit that best explains the sequence structure. The algorithm tries all
possible repeat unit sizes and selects the best segmentation based on:

- Highest 3' flanking repeat count (most important)

- Highest internal repeat count

- Lowest spacer length (prefer clean repeats)

- Lowest unit length (prefer simpler units)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple AT repeat
result <- seg_simple("d", "ATATAT", "ATATGG")
print(result)

# CG repeat
result <- seg_simple("d", "CGCGCG", "CGCGAA")
print(result)
} # }
```
