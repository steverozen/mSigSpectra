# Check and, if possible, correct the chromosome names in a VCF data.frame

Harmonizes the `CHROM` column of `vcf.df` with the chromosome names used
by `ref.genome` (a BSgenome object). Adds or strips the `chr` prefix as
needed; for human and mouse, translates `23/24` or `20/21` to `X/Y`.

## Usage

``` r
check_and_fix_chrom_names(vcf.df, ref.genome, name.of.VCF = NULL)
```

## Arguments

- vcf.df:

  A VCF as a data frame with a `CHROM` column.

- ref.genome:

  A BSgenome object (e.g. `BSgenome.Hsapiens.UCSC.hg38`).

- name.of.VCF:

  Name of the VCF file (for warning/error messages).

## Value

A character vector of chromosome names that can be used as a replacement
for `vcf.df$CHROM`. Errors if reconciliation is impossible.
