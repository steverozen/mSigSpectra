# Annotate a VCF with flanking sequence context and transcript strand

`annotate_vcf()` is the high-level annotator used by
[`vcf_to_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_catalog.md).
It dispatches on `variant_type`:

## Usage

``` r
annotate_vcf(
  vcf,
  ref_genome,
  variant_type = c("auto", "SBS", "DBS", "ID"),
  trans_ranges = NULL,
  seq_context_width = 10L,
  name_of_vcf = NULL
)
```

## Arguments

- vcf:

  A VCF as a data.frame / data.table with `CHROM`, `POS`, `REF`, `ALT`
  columns.

- ref_genome:

  A BSgenome object or a character alias accepted by
  [`normalize_genome_arg()`](https://steverozen.github.io/mSigSpectra/reference/normalize_genome_arg.md).

- variant_type:

  One of `"auto"`, `"SBS"`, `"DBS"`, `"ID"`. For `"auto"`, rows are
  routed by REF/ALT length.

- trans_ranges:

  Optional transcript-ranges `data.table`. If `NULL`, the shipped table
  is used when available.

- seq_context_width:

  Width (per side) of the flanking-sequence window; default 10 (produces
  a 21-base window).

- name_of_vcf:

  Optional VCF name used in warnings.

## Value

A `data.table` with the annotation columns added.

## Details

- `"SBS"`: add flanking sequence context
  ([`add_seq_context()`](https://steverozen.github.io/mSigSpectra/reference/add_seq_context.md))
  and transcript strand
  ([`add_transcript_strand()`](https://steverozen.github.io/mSigSpectra/reference/add_transcript_strand.md)
  if `trans_ranges` is supplied or shipped for `ref_genome`).

- `"DBS"`: same as SBS; downstream
  [`vcf_to_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_catalog.md)
  will require the tetranucleotide context to be free of `N`.

- `"ID"`: placeholder — indel annotation will require indel
  justification and repeat-context detection. Not yet implemented in
  this port.

- `"auto"`: classify each row by REF/ALT length and annotate each subset
  appropriately. Only the SBS / DBS paths are wired up for now.
