# Remove rows that look like stray "#CHROM" header repeats

Occasionally a VCF body contains a spurious row whose `CHROM` column
equals `"#CHROM"` (e.g. when a concatenated VCF keeps the header line of
each input). Drop those rows and record them in `discarded.variants`.

## Usage

``` r
remove_rows_with_pound_sign(df, name_of_vcf = NULL)
```
