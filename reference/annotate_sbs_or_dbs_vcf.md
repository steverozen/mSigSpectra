# Annotate an SBS or DBS VCF with flanking sequence context and transcript strand

Adds flanking sequence context via
[`add_seq_context()`](https://steverozen.github.io/mSigSpectra/reference/add_seq_context.md)
and, when a transcript-ranges table is available, transcript strand via
[`add_transcript_strand()`](https://steverozen.github.io/mSigSpectra/reference/add_transcript_strand.md).
The same pipeline is appropriate for SBS and DBS; DBS-specific
validation (e.g. no `N` in the tetranucleotide context) happens
downstream in
[`vcf_to_dbs_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_dbs_catalog.md).

## Usage

``` r
annotate_sbs_or_dbs_vcf(
  vcf,
  ref_genome,
  trans_ranges = NULL,
  seq_context_width = 10L,
  name_of_vcf = NULL
)
```

## Arguments

- vcf:

  A VCF as a data.frame / data.table with `CHROM`, `POS`, `REF`, `ALT`
  columns. For SBS, `REF` and `ALT` are single bases; for DBS, both are
  two bases.

- ref_genome:

  A BSgenome object or a character alias accepted by
  [`normalize_genome_arg()`](https://steverozen.github.io/mSigSpectra/reference/normalize_genome_arg.md).

- trans_ranges:

  Optional transcript-ranges `data.table`. If `NULL`, the shipped table
  is used when available for `ref_genome`.

- seq_context_width:

  Width (per side) of the flanking-sequence window (default 10 → 21-base
  window).

- name_of_vcf:

  Optional VCF name used in warnings.

## Value

A list with:

- `annotated.vcf`: the input VCF with new columns `seq.<N>bases` and
  (when transcript ranges are available) `trans.strand` / `bothstrand`.

- `discarded.variants`: currently always `NULL` for this path (included
  for symmetry with
  [`annotate_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_id_vcf.md)).
