# Add flanking sequence context to a VCF data frame

Extracts `seq_context_width` bases upstream and downstream of each
variant position from the reference genome and attaches them as a new
column named `seq.<N>bases` where `N = 2 * seq_context_width + 1`.

## Usage

``` r
add_seq_context(df, ref_genome, seq_context_width = 10, name_of_vcf = NULL)
```

## Arguments

- df:

  A VCF as a data frame / data.table with columns `CHROM` and `POS`.

- ref_genome:

  A BSgenome object or a character identifier accepted by
  [`normalize_genome_arg()`](https://steverozen.github.io/mSigSpectra/reference/normalize_genome_arg.md).

- seq_context_width:

  Number of flanking bases on each side (default 10, producing a 21-base
  window).

- name_of_vcf:

  Optional VCF name used in warning/error messages.

## Value

`df` with an added character column `seq.<N>bases`.
