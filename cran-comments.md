# Resubmission

This is a resubmission of mSigSpectra (now 0.1.3) addressing the CRAN
review of 0.1.2 (Leonore Hochhauser, 2026-09-15). Each point and its
resolution:

* **Explain all acronyms in the Description.** VCF, SBS, DBS, ID
  (indel), and PDF are now spelled out on first use, in both Title and
  Description.
* **Quote package / software names.** `'ICAMS'`, `'shiny'`, and
  `'mSigPlot'` are single-quoted in the Description.
* **License.** `License: GPL-3`. The `+ file LICENSE` component and the
  `LICENSE` file have been removed.
* **Missing `\value` tags** in `check_and_remove_discarded_variants.Rd`,
  `is_catalog.Rd`, and `subset_catalog.Rd`. Added, describing the class
  and meaning of each return value. `\value` was also added to every
  internal helper's Rd file.
* **`print()` in `R/quick_check_vcf.R`.** The multiple-ALT warning text
  is now built with `sprintf()` / `paste()` and passed to `warning()`.
  There are no remaining `print()` / `cat()` calls in the package code.

mSigSpectra reads variant call files (VCFs) in a caller-agnostic way,
annotates variants with flanking sequence context and transcriptional
strand, and builds mutational-spectrum catalogs (SBS, DBS, indel) at
several resolutions (SBS96/192/1536, DBS78/136/144, ID83/89/166/476)
for GRCh37, GRCh38, and GRCm38.

## R CMD check results

0 ERRORs, 0 WARNINGs across all five environments in the
r-lib/actions/check-r-package GitHub Actions matrix
(macOS-latest release, Windows-latest release, Ubuntu-latest devel /
release / oldrel-1; run
https://github.com/steverozen/mSigSpectra/actions/runs/35125181543),
and 0 ERRORs / 0 WARNINGs on a local `R CMD check --as-cran` run on
the built tarball with the four CRAN incoming-feasibility env vars
enabled (`_R_CHECK_CRAN_INCOMING_`,
`_R_CHECK_CRAN_INCOMING_REMOTE_`,
`_R_CHECK_CRAN_INCOMING_CHECK_FILE_URIS_`,
`_R_CHECK_CRAN_INCOMING_USE_ASPELL_`).

## Notes seen on CI

Each of the five CI jobs reported a single NOTE from the CRAN
incoming-feasibility check:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Steve Rozen <steverozen@pm.me>'
New submission
No suitable spell-checker program found
```

The "New submission" line is expected, this is the first CRAN
submission of mSigSpectra. The "No suitable spell-checker program
found" line is emitted because the GH Actions runners do not have
aspell installed; it is not a defect of the package.

## Notes seen on local `R CMD check --as-cran`

The local check (with aspell installed) additionally reports:

```
Possibly misspelled words in DESCRIPTION:
  DBS (31:54, 33:22, 33:29, 33:37)
  Rozen (37:58)
  SBS (31:20, 32:67, 33:5, 33:13)
  VCF (28:41)
  al (37:67)
  et (37:64)
  indels (32:31)
  transcriptional (29:59)
```

All of these are intentional. `VCF`, `SBS`, `DBS`, and `indels` are
now defined on first use in the Description and then used as
abbreviations. `transcriptional` is a standard term. `Rozen` is the
maintainer's surname and `et al` is part of the citation.

Two further local NOTEs are artifacts of the local machine, not the
package: "Compilation used the following non-portable flag(s):
'-march=native'" comes from the local `~/.R/Makevars`, and "no command
'tidy' found" reflects a missing HTML Tidy binary.

## Reverse dependencies

There are no reverse dependencies on CRAN (this is a new package).
