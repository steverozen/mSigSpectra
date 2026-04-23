# Split a mixed-mutation VCF into SBS / DBS / ID sub-tables

Classifies each row by `REF`/`ALT` length alone:

- SBS (single base substitution): `nchar(REF) == 1 && nchar(ALT) == 1`

- DBS (double base substitution): `nchar(REF) == 2 && nchar(ALT) == 2`

- ID (indel): `nchar(REF) != nchar(ALT)`

- Other (e.g. 3+bp substitutions): discarded with a recorded reason.

## Usage

``` r
split_vcf(vcf, name_of_vcf = NULL)
```

## Arguments

- vcf:

  A VCF as a data.frame / data.table with at least `REF` and `ALT`
  columns.

- name_of_vcf:

  Optional VCF name used in warning / error messages.

## Value

A list with elements

- `SBS`: `data.table` of SBS rows.

- `DBS`: `data.table` of DBS rows.

- `ID`: `data.table` of indel rows.

- `discarded`: `data.table` of rows that did not fit any of the above,
  with a `discarded.reason` column — **`NULL` when no rows were
  discarded**.

## Details

This is **caller-agnostic**. In particular, mSigSpectra does *not* merge
adjacent SBSs into DBSs based on VAF similarity (which ICAMS did via
`SplitOneVCF` for Strelka-style VCFs). Users who want that behavior
should apply it as a post-processing step with their own VAF column.
