# mSigSpectra 0.1.3

* Changes requested in the CRAN review of 0.1.2:
  * `License` field is now `GPL-3` and the redundant `LICENSE` file
    has been removed.
  * Acronyms (VCF, SBS, DBS, ID, PDF) are spelled out in DESCRIPTION
    and software names are single-quoted.
  * Added `\value` sections to the documentation of
    `check_and_remove_discarded_variants()`, `is_catalog()`, and
    `subset_catalog()`, and to all internal helper documentation.
  * `quick_check_vcf()` no longer uses `print()` to format the
    multiple-ALT warning.

# mSigSpectra 0.1.2

* Initial CRAN submission.
* Provides a four-step pipeline for building mutational-spectrum
  catalogs from VCF files: `read_vcf()`, `split_vcf()`,
  `annotate_sbs_or_dbs_vcf()` / `annotate_id_vcf()`, and
  `vcf_to_sbs_catalog()` / `vcf_to_dbs_catalog()` /
  `vcf_to_id_catalog()`.
* Supports SBS, DBS, and indel (ID) catalogs for GRCh37, GRCh38,
  and GRCm38 reference genomes or any BSgenome genome.
  Adds transcript-strand annotation where transcript-ranges info is available.
* Ships example Strelka and Mutect VCFs in `inst/extdata/` and a
  vignette demonstrating the full pipeline.

