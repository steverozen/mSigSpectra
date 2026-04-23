# Package index

## All functions

- [`add_seq_context()`](https://steverozen.github.io/mSigSpectra/reference/add_seq_context.md)
  : Add flanking sequence context to a VCF data frame
- [`add_transcript_strand()`](https://steverozen.github.io/mSigSpectra/reference/add_transcript_strand.md)
  : Annotate a VCF data frame with transcript strand information
- [`all.abundance`](https://steverozen.github.io/mSigSpectra/reference/all.abundance.md)
  : K-mer abundances for density calculations
- [`annot_vcf_to_476_catalog()`](https://steverozen.github.io/mSigSpectra/reference/annot_vcf_to_476_catalog.md)
  : Convert an annotated indel VCF to a Koh 476-category catalog
- [`annot_vcf_to_83_catalog()`](https://steverozen.github.io/mSigSpectra/reference/annot_vcf_to_83_catalog.md)
  : Convert an annotated indel VCF to a COSMIC 83-category catalog
- [`annot_vcf_to_89_catalog()`](https://steverozen.github.io/mSigSpectra/reference/annot_vcf_to_89_catalog.md)
  : Convert an annotated indel VCF to a Koh 89-category catalog
- [`annotate_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_id_vcf.md)
  : Annotate an in-memory ID (indel) VCF with sequence context,
  transcript strand, and COSMIC / Koh indel categories
- [`annotate_sbs_or_dbs_vcf()`](https://steverozen.github.io/mSigSpectra/reference/annotate_sbs_or_dbs_vcf.md)
  : Annotate an SBS or DBS VCF with flanking sequence context and
  transcript strand
- [`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)
  : Turn a numeric matrix into a mutational-spectrum catalog
- [`catalog.row.order`](https://steverozen.github.io/mSigSpectra/reference/catalog.row.order.md)
  : Canonical row orders for catalogs
- [`catalog_attrs()`](https://steverozen.github.io/mSigSpectra/reference/catalog_attrs.md)
  : Report the attributes of an mSigSpectra catalog
- [`categorize_1_justified_indel()`](https://steverozen.github.io/mSigSpectra/reference/categorize_1_justified_indel.md)
  : Given a indel and its sequence context, categorize it
- [`cbind_catalogs()`](https://steverozen.github.io/mSigSpectra/reference/cbind_catalogs.md)
  : Combine catalogs across samples (column-bind)
- [`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md)
  : Check a VCF for common variant-level problems and remove the
  offenders
- [`collapse_catalog()`](https://steverozen.github.io/mSigSpectra/reference/collapse_catalog.md)
  : Collapse a higher-resolution catalog to a lower-resolution one
- [`infer_trans_ranges()`](https://steverozen.github.io/mSigSpectra/reference/infer_trans_ranges.md)
  : Infer transcript ranges for a reference genome
- [`is_catalog()`](https://steverozen.github.io/mSigSpectra/reference/is_catalog.md)
  : Check whether an object looks like an mSigSpectra catalog
- [`justify_id_vcf()`](https://steverozen.github.io/mSigSpectra/reference/justify_id_vcf.md)
  : Add sequence context and transcript information to an in-memory ID
  (insertion/deletion) VCF, and confirm that they match the given
  reference genome
- [`justify_indel()`](https://steverozen.github.io/mSigSpectra/reference/justify_indel.md)
  : Move the notional position of a deletion or insertion as far left as
  possible.
- [`read_catalog()`](https://steverozen.github.io/mSigSpectra/reference/read_catalog.md)
  : Read a mutational-spectrum catalog from a file
- [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md)
  : Read a VCF file into a data.table, caller-agnostically
- [`read_vcfs()`](https://steverozen.github.io/mSigSpectra/reference/read_vcfs.md)
  : Read multiple VCF files
- [`revc()`](https://steverozen.github.io/mSigSpectra/reference/revc.md)
  : Reverse complement a character vector of DNA strings
- [`seg_simple()`](https://steverozen.github.io/mSigSpectra/reference/seg_simple.md)
  : Segment a single indel sequence using Rcpp interface
- [`segment_simple_cpp()`](https://steverozen.github.io/mSigSpectra/reference/segment_simple_cpp.md)
  : Segment a single indel using Rcpp interface
- [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)
  : Split a mixed-mutation VCF into SBS / DBS / ID sub-tables
- [`subset_catalog()`](https://steverozen.github.io/mSigSpectra/reference/subset_catalog.md)
  : Subset a catalog while preserving attributes
- [`trans.ranges.GRCh37`](https://steverozen.github.io/mSigSpectra/reference/trans.ranges.md)
  [`trans.ranges.GRCh38`](https://steverozen.github.io/mSigSpectra/reference/trans.ranges.md)
  [`trans.ranges.GRCm38`](https://steverozen.github.io/mSigSpectra/reference/trans.ranges.md)
  : Transcript ranges for transcriptional strand annotation
- [`transform_catalog()`](https://steverozen.github.io/mSigSpectra/reference/transform_catalog.md)
  : Transform a catalog between counts and density
- [`vcf_to_dbs_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_dbs_catalog.md)
  : Build a DBS mutational-spectrum catalog from an annotated DBS VCF
- [`vcf_to_id_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_id_catalog.md)
  : Build an ID (indel) mutational-spectrum catalog from an annotated ID
  VCF
- [`vcf_to_sbs_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_sbs_catalog.md)
  : Build an SBS mutational-spectrum catalog from an annotated SBS VCF
- [`write_catalog()`](https://steverozen.github.io/mSigSpectra/reference/write_catalog.md)
  : Write a mutational-spectrum catalog to a file
