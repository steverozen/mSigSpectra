# mSigSpectra

Build mutational-spectrum catalogs from VCF files.

mSigSpectra is the **lean, caller-agnostic successor to
[ICAMS](https://github.com/steverozen/ICAMS)** — it keeps the numerical
core (VCF reading, annotation, indel justification, and SBS / DBS / ID
catalog construction) and drops everything else. **Plotting is now
handled by the separate
[`mSigPlot`](https://github.com/steverozen/mSigPlot) package.** Shiny,
PDF / zip reporting, and per-caller VAF extraction are intentionally out
of scope.

## What it does

    read_vcf()  →  split_vcf()  →  annotate_vcf()  →  vcf_to_catalog()

- [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md)
  — reads any VCF file body. **Caller-agnostic:** it does not parse
  `FORMAT` / sample columns or extract VAF.
- [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)
  — partitions a VCF into SBS / DBS / ID sub-tables by `REF`/`ALT`
  length alone.
- `annotate_vcf()` — adds flanking sequence context
  ([`BSgenome::getSeq`](https://rdrr.io/pkg/Biostrings/man/getSeq.html)),
  transcript strand
  ([`GenomicRanges::findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  against shipped GENCODE tables), and — for indels — left-justifies and
  categorizes them (COSMIC-83 / Koh-89 / Koh-476).
- [`vcf_to_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_catalog.md)
  — produces a single catalog of any of these types: `SBS96`, `SBS192`,
  `SBS1536`, `DBS78`, `DBS136`, `DBS144`, `ID83`, `ID89`, `ID166`,
  `ID476`.

The output is a plain numeric matrix with attributes (`type`,
`counts_or_density`, `ref_genome`, `region`, `abundance`). **No S3
classes.** Functions key off `attr(x, "type")` via explicit checks.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("steverozen/mSigSpectra")
```

You also need the BSgenome package(s) for your reference genome(s):

``` r
BiocManager::install(c(
  "BSgenome.Hsapiens.1000genomes.hs37d5",   # GRCh37 / hg19
  "BSgenome.Hsapiens.UCSC.hg38",            # GRCh38 / hg38
  "BSgenome.Mmusculus.UCSC.mm10"            # GRCm38 / mm10
))
```

For plotting, install mSigPlot:

``` r
remotes::install_github("steverozen/mSigPlot")
```

## Quick example

``` r
library(mSigSpectra)
library(mSigPlot)

vcf <- read_vcf("my_sample.vcf", filter = "PASS")
parts <- split_vcf(vcf)

ann_sbs <- annotate_vcf(parts$SBS, ref_genome = "GRCh38",
                        variant_type = "SBS")
sbs96 <- vcf_to_catalog(ann_sbs, type = "SBS96",
                        ref_genome = "GRCh38", region = "genome",
                        sample_name = "my_sample")

plot_SBS96(sbs96)
```

See
[`vignettes/mSigSpectra.Rmd`](https://steverozen.github.io/mSigSpectra/vignettes/mSigSpectra.Rmd)
for a full walk-through that reads a real VCF, generates SBS / DBS / ID
catalogs, writes them to disk, and plots each one with mSigPlot.

## Migrating from ICAMS

| ICAMS                                                       | mSigSpectra                                                                                                                                                                       |
|-------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ReadVCFs()` / `SimpleReadVCF()`                            | [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md), [`read_vcfs()`](https://steverozen.github.io/mSigSpectra/reference/read_vcfs.md)                  |
| `SplitListOfVCFs()`                                         | [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)                                                                                                  |
| `AnnotateSBSVCF()` / `AnnotateDBSVCF()` / `AnnotateIDVCF()` | `annotate_vcf(variant_type = ...)`                                                                                                                                                |
| `VCFsToCatalogs()`                                          | `vcfs_to_catalogs()`                                                                                                                                                              |
| `as.catalog()`                                              | [`as_catalog()`](https://steverozen.github.io/mSigSpectra/reference/as_catalog.md)                                                                                                |
| `TransformCatalog()`                                        | [`transform_catalog()`](https://steverozen.github.io/mSigSpectra/reference/transform_catalog.md)                                                                                  |
| `ReadCatalog()` / `WriteCatalog()`                          | [`read_catalog()`](https://steverozen.github.io/mSigSpectra/reference/read_catalog.md) / [`write_catalog()`](https://steverozen.github.io/mSigSpectra/reference/write_catalog.md) |
| `PlotCatalog()`                                             | (now in [`mSigPlot::plot_SBS96()`](https://steverozen.github.io/mSigPlot/reference/bar_plots.html) etc.)                                                                          |

UpperCamelCase is gone. S3 catalog classes (`SBS96Catalog`,
`IndelCatalog`, …) are gone — catalogs are plain matrices with
attributes. The `variant.caller` argument and the per-caller VAF
extractors (`GetStrelkaVAF`, `GetMutectVAF`, etc.) are gone.

## Gotchas you have to handle yourself

mSigSpectra deliberately stays out of several judgment calls that ICAMS
made for you. **The reader returns the rows the VCF actually contains.**
You are responsible for sanitizing them before downstream analysis if
the caller / VCF source has any of the following quirks:

### 1. Adjacent SBSs are *not* merged into DBSs

Some callers (notably Strelka) emit a real DBS as two adjacent SBSs at
positions *p* and *p+1*. ICAMS’s Strelka-specific code path tried to
detect these by VAF similarity and merged them into a single DBS row.
**mSigSpectra never does this** — the result is one fewer DBS and two
extra SBSs per pair. If your VCFs come from such a caller and you need
DBSs, merge adjacent SBSs yourself (using your own VAF column) **before
calling
[`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)**.

### 2. Multi-allelic sites: `ALT = "G,T"` style

Tri-allelic and tetra-allelic sites can be encoded as a single VCF row
with multiple comma-separated alleles in `ALT` (`A → G,T`), or as
multiple rows at the same `(CHROM, POS)`, or with the
`bcftools norm -m -` “split” representation. mSigSpectra:

- **Multi-allelic-on-one-row** rows (`ALT` contains `,`) —
  [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md)
  passes them through.
  [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)
  will not classify them as SBS / DBS / ID since `nchar(ALT)` is 3+.
  They land silently in `split_vcf()$discarded`. **If you want to
  include them, split each alternate allele into its own row first**
  (e.g. with `bcftools norm -m -any`).
- **Same-position multi-row** rows are kept by
  [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md).
  The optional
  [`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md)
  will discard them with reason “Variant with same CHROM, POS, REF but
  different ALT” — call it explicitly if you want this behavior.

### 3. “Complex” indels and equivalent-length substitutions

A complex change like `REF = "AG"`, `ALT = "TC"` (two-base substitution
where neither base matches) can equally well be represented as two
adjacent SBSs (`A → T` and `G → C`), or as a single DBS, or as one
deletion plus one insertion. There is **no community standard.**
mSigSpectra:

- **`REF[1] != ALT[1]` and lengths differ** — discarded by
  [`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md)
  as “Complex indel”; left in place by the bare
  [`read_vcf()`](https://steverozen.github.io/mSigSpectra/reference/read_vcf.md) +
  [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)
  path (where `nchar(REF) != nchar(ALT)` flags it as ID even though it
  does not conform to the ICAMS / COSMIC indel convention of a shared
  first base).
- **3-bp+ same-length substitutions** like `ACT → TGA` — discarded by
  [`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md)
  as “Variant involves three or more nucleotides”; classified into
  `split_vcf()$discarded` by the
  [`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)
  path. If you have a curated convention for these, apply it before
  catalog construction.

### 4. Invalid DBSs sharing a base

Some callers emit `REF = "TA"`, `ALT = "TT"` (the `T` at position 1 is
unchanged). These are degenerate DBSs — really single-base substitutions
in disguise.
[`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md)
discards them as “Wrong DBS variant”;
[`split_vcf()`](https://steverozen.github.io/mSigSpectra/reference/split_vcf.md)
will classify them as DBS unless you filter first.

### 5. Ambiguous reference bases

Variants with non-`{A,C,G,T}` REF bases (e.g. `N`) cannot be
canonicalized.
[`check_and_remove_discarded_variants()`](https://steverozen.github.io/mSigSpectra/reference/check_and_remove_discarded_variants.md)
removes them;
[`vcf_to_catalog()`](https://steverozen.github.io/mSigSpectra/reference/vcf_to_catalog.md)
ultimately drops them when the pentanucleotide context contains `N` and
warns.

### 6. FILTER conventions vary by caller

The default `read_vcf(filter = TRUE)` keeps rows where `FILTER` is in
`{"PASS", ".", ""}` — the union of common caller defaults. **This is not
bit-exact compatible with any single caller** — Strelka and Mutect keep
only `"PASS"`, Freebayes uses `"."`. If you need parity with one
specific caller, pass `filter = "PASS"` (or your caller’s exact status
string) explicitly.

### Suggested defensive pipeline

``` r
vcf <- read_vcf(file, filter = "PASS")
clean <- check_and_remove_discarded_variants(vcf, name_of_vcf = file)
parts <- split_vcf(clean$df, name_of_vcf = file)
# Inspect clean$discarded.variants and parts$discarded for anything you
# want to recover or re-encode before the annotate / catalog step.
```

## Status

- Version 0.0.1.
- `R CMD check`: 0 errors, 0 warnings.
- 235 testthat assertions across 13 test files.
- SBS96 catalog construction is bit-exact against ICAMS for Mutect VCFs;
  diverges from ICAMS’s Strelka path by the documented
  no-adjacency-merge design choice.
- DBS catalog construction passes structural / shape tests; full numeric
  parity vs ICAMS not yet verified (planned).
- ID83 / ID89 / ID476 catalogs are wired through end-to-end on real
  Strelka ID VCFs.
- [`transform_catalog()`](https://steverozen.github.io/mSigSpectra/reference/transform_catalog.md)
  (counts ↔︎ density) and
  [`collapse_catalog()`](https://steverozen.github.io/mSigSpectra/reference/collapse_catalog.md)
  (SBS1536 → SBS96, SBS192 → SBS96, DBS144 → DBS78) implemented.
- [`read_catalog()`](https://steverozen.github.io/mSigSpectra/reference/read_catalog.md)
  /
  [`write_catalog()`](https://steverozen.github.io/mSigSpectra/reference/write_catalog.md)
  round-trip ICAMS-native CSV losslessly; SigProfiler input formats
  parse for SBS96, SBS1536, ID83.

## License

GPL-3.

## Citation

If you use mSigSpectra in published work, please cite the original ICAMS
papers:

- Boot et al., *Genome Research* 2018,
  [doi:10.1101/gr.230219.117](https://doi.org/10.1101/gr.230219.117).
- Boot et al., *Genome Research* 2020,
  [doi:10.1101/gr.255620.119](https://doi.org/10.1101/gr.255620.119).
