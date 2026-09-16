# Changelog

## mSigSpectra 0.1.3

- Changes requested in the CRAN review of 0.1.2:
  - `License` field is now `GPL-3` and the redundant `LICENSE` file has
    been removed.
  - Acronyms (VCF, SBS, DBS, ID, PDF) are spelled out in DESCRIPTION and
    software names are single-quoted.
  - Added `\value` sections to the documentation of
    [`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md),
    [`is_catalog()`](https://steverozen.github.io/mSigSpectra/reference/is_catalog.md),
    and
    [`subset_catalog()`](https://steverozen.github.io/mSigSpectra/reference/subset_catalog.md),
    and to all internal helper documentation.
  - [`quick_check_vcf()`](https://steverozen.github.io/mSigSpectra/reference/quick_check_vcf.md)
    no longer uses [`print()`](https://rdrr.io/r/base/print.html) to
    format the multiple-ALT warning.

## mSigSpectra 0.1.2

- Initial CRAN submission.
- Provides a four-step pipeline for building mutational-spectrum
  catalogs from VCF files:
  [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md),
  [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md),
  [`annotate_sbs_or_dbs_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_sbs_or_dbs_vcf.md)
  /
  [`annotate_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_id_vcf.md),
  and
  [`vcf_to_sbs_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_sbs_catalog.md)
  /
  [`vcf_to_dbs_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_dbs_catalog.md)
  /
  [`vcf_to_id_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_id_catalog.md).
- Supports SBS, DBS, and indel (ID) catalogs for GRCh37, GRCh38, and
  GRCm38 reference genomes or any BSgenome genome. Adds
  transcript-strand annotation where transcript-ranges info is
  available.
- Ships example Strelka and Mutect VCFs in `inst/extdata/` and a
  vignette demonstrating the full pipeline.
