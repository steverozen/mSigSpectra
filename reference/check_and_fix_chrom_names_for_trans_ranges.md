# Check and, if possible, correct the chromosome names in a trans.ranges table

Harmonizes `trans.ranges$chrom` with the chromosome naming used by
`vcf.df$CHROM`, adding or stripping `chr` as needed and translating
organism-specific numeric X/Y encodings.

## Usage

``` r
check_and_fix_chrom_names_for_trans_ranges(
  trans.ranges,
  vcf.df,
  ref.genome,
  name.of.VCF = NULL
)
```

## Arguments

- trans.ranges:

  A `data.table` of transcript ranges (see
  [trans.ranges](https://steverozen.github.io/mSigSpectra/reference/trans.ranges.md)).

- vcf.df:

  A VCF as a data frame with a `CHROM` column.

- ref.genome:

  A BSgenome object used to determine organism.

- name.of.VCF:

  Name of the VCF file.

## Value

A character vector of chromosome names that can be used as a replacement
for `trans.ranges$chrom`. Errors if reconciliation is impossible.
