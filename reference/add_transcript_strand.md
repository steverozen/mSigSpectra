# Annotate a VCF data frame with transcript strand information

For each variant, finds overlapping transcripts in `trans_ranges` via
[`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/GenomicRanges/man/findOverlaps-methods.html)
with `type = "within"`, and appends columns `trans.start.pos`,
`trans.end.pos`, `trans.strand`, `trans.Ensembl.gene.ID`,
`trans.gene.symbol`, plus `bothstrand` (TRUE if the variant falls on
transcripts from both strands) and `count` (number of overlapping
transcripts).

## Usage

``` r
add_transcript_strand(df, ref_genome, trans_ranges = NULL, name_of_vcf = NULL)
```

## Arguments

- df:

  A VCF as a data frame / data.table with columns `CHROM`, `POS`, `ALT`.

- ref_genome:

  A BSgenome object or a character identifier accepted by
  [`normalize_genome_arg()`](https://steverozen.github.io/mSigSpectra/reference/normalize_genome_arg.md).

- trans_ranges:

  Optional `data.table` of transcript ranges with columns `chrom`,
  `start`, `end`, `strand`, `Ensembl.gene.ID`, `gene.symbol`. If `NULL`,
  the shipped table for `ref_genome` is used.

- name_of_vcf:

  Optional VCF name used in warning/error messages.

## Value

A `data.table` with the annotation columns added (left-join semantics:
variants outside any transcript have `NA` values).

## Details

If `trans_ranges` is `NULL` and no shipped table is available for the
given `ref_genome`, returns `df` as a `data.table` unchanged (no strand
columns added). If `df` has zero rows, it is returned unchanged.
