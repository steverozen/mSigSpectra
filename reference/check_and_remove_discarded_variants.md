# Check a VCF for common variant-level problems and remove the offenders

Removes:

- Rows with identical REF and ALT.

- Stray `#CHROM` header repeats.

- Duplicated (CHROM, POS, REF, ALT) rows (keeping one copy).

- Multiple-ALT rows (comma-separated ALT field).

- Non-standard chromosome names (or those outside `chr_names_to_process`
  when supplied).

- Substitutions of length \> 2 (e.g. `ACT>TGA`).

- Complex indels (`REF[1] != ALT[1]`).

- Wrong DBS rows where REF and ALT share a base at the same position.

- Variants whose REF base is not in `{A, C, G, T}`.

## Usage

``` r
check_and_remove_discarded_variants(
  vcf,
  name_of_vcf = NULL,
  chr_names_to_process = NULL
)
```

## Arguments

- vcf:

  A VCF as a data.frame / data.table.

- name_of_vcf:

  Optional name, used in warning messages.

- chr_names_to_process:

  Optional character vector of chromosome names to keep (overrides the
  default non-standard-contig filter).

## Details

Each discarded row gains a `discarded.reason` column. Returns a list
with `df` (retained rows) and optionally `discarded.variants` (discarded
rows).
