# Justify indels in a VCF and update positions accordingly

For each indel in the input VCF, this function:

1.  Justifies the indel position (moves it as far left as possible)

2.  Calculates how much the position was shifted

3.  Updates the POS column by decrementing it by the shift amount

4.  Updates the seq.context.width column accordingly

5.  Adds a new column 'pos_shift' showing how much each position was
    moved

6.  Categorizes the justified indel

## Usage

``` r
justify_indels_in_id_vcf_with_contexts(vcf, explain_indels = 1)
```

## Arguments

- vcf:

  A data.frame representing a VCF, which must contain the following
  columns:

  - seq.context: Ample surrounding sequence on each side of the variants

  - REF: The reference alleles (includes one unaltered base at the
    start)

  - ALT: The alternative alleles

  - POS: The genomic positions

  - seq.context.width: The width of seq.context to the left

- explain_indels:

  0 stay silent, if 1, print explanation of there was a change in POS,
  if 3, aways print

## Value

A data.frame with:

- All original columns from the input vcf

- Updated POS values (decremented by pos_shift)

- Updated REF and ALT values if the positon of the indel has been moved

- Updated seq.context.width values (decremented by pos_shift)

- New column 'pos_shift': Amount the position was moved (usually 0)

## Important invariant

The seq.context string is never modified. Both POS and seq.context.width
are decremented by the same amount (pos_shift), which maintains the
relationship: genomic position POS corresponds to position
(seq.context.width + 1) in seq.context.

For example, if seq.context was extracted from genomic positions (POS -
seq.context.width, POS + var.width + seq.context.width), then after
justification by shift S:

- seq.context remains unchanged

- POS becomes POS - S

- seq.context.width becomes seq.context.width - S

- The variant is now at position (seq.context.width - S) + 1 in
  seq.context
