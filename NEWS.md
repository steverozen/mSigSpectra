# mSigSpectra 0.1.1

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
