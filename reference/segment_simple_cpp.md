# Segment a single indel using Rcpp interface

This function provides direct C++ interface for indel segmentation
without calling an external binary via system2.

## Usage

``` r
segment_simple_cpp(ins_or_del, string, context)
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
