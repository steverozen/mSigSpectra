# Remove variants with duplicated CHROM+POS

Two passes:

1.  Rows with identical (CHROM, POS, REF, ALT): keep one copy, discard
    the rest.

2.  Rows sharing (CHROM, POS, REF) but different ALT: discard all
    (treated as unresolved multiallelic / inconsistent records).

## Usage

``` r
remove_rows_with_duplicated_chrom_and_pos(df, name_of_vcf = NULL)
```
