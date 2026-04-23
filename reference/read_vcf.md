# Read a VCF file into a data.table, caller-agnostically

Reads the *body* of a VCF file (lines after `#CHROM`) into a
`data.table`. The resulting table has whatever columns the VCF has
(`CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`, and
optionally `FORMAT` plus one or more sample columns), with `#CHROM`
renamed to `CHROM`.

## Usage

``` r
read_vcf(
  file,
  filter = TRUE,
  backend = c("fread", "vcfppR"),
  name_of_vcf = NULL
)
```

## Arguments

- file:

  Path or URL to the VCF file.

- filter:

  Controls which rows are kept based on the `FILTER` column:

  - `TRUE` (default): keep rows where `FILTER %in% c("PASS", ".", "")` —
    the union of passing values across common callers.

  - `FALSE` or `NULL`: keep all rows.

  - A character vector: keep rows where `FILTER %in% filter`. Rows with
    no `FILTER` column are always kept; a warning is emitted when
    `filter` is non-trivial in that case.

- backend:

  VCF reader backend:

  - `"fread"` (default):
    [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).
    Handles uncompressed and gzipped files; does not handle
    bgzipped/tabix.

  - `"vcfppR"`:
    [`vcfppR::vcftable()`](https://rdrr.io/pkg/vcfppR/man/vcftable.html)
    — handles bgzipped, tabix, and remote files; splits multi-allelic
    records. The `vcfppR` package is in `Suggests`; an informative error
    is raised if it is not installed.

- name_of_vcf:

  Optional name for the VCF, used only for warning / error messages.
  Defaults to the filename with extension stripped.

## Value

A `data.table` with one row per variant.

## Details

**Caller-agnostic.** `read_vcf()` does not know or care which variant
caller produced the VCF. It does not parse `FORMAT`/sample columns and
does not extract VAF or read depth. The only caller-dependent semantics
is the default value of the `filter` argument (see below).
