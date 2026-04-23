# Read a mutational-spectrum catalog from a file

Auto-detects the file format from its first row. Recognized formats:

- **ICAMS native** CSV — first N columns are header columns (mutation
  type, context, etc), remaining columns are per-sample counts /
  densities.

- **SigProfiler** TSV / CSV — single header column with bracketed
  mutation labels like `A[C>A]A`, `AA[C>A]AA`, or PCAWG indel codes like
  `1:Del:C:0`.

- **COSMIC** CSV — SBS96 uses the ICAMS-external row-header format;
  SBS192 uses a stranded variant of the same; DBS78 / ID83 also share
  their SigProfiler-style layout.

## Usage

``` r
read_catalog(
  file,
  ref_genome = NULL,
  region = "unknown",
  counts_or_density = "counts",
  format = c("auto", "ICAMS", "SigProfiler", "COSMIC")
)
```

## Arguments

- file:

  Path to the catalog file.

- ref_genome:

  Optional BSgenome object or alias; stored as attribute.

- region:

  One of `"genome"`, `"exome"`, `"transcript"`, `"unknown"`.

- counts_or_density:

  One of `"counts"`, `"density"`, `"counts.signature"`,
  `"density.signature"`.

- format:

  `"auto"` (default), `"ICAMS"`, `"SigProfiler"`, or `"COSMIC"`. In
  `"auto"` mode the format is inferred from the file's first data row.

## Value

A catalog matrix with attributes (see
[`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)).

## Details

The matrix is reordered to the canonical rownames for its type and
wrapped in a catalog via
[`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md).
