# Submission

This is the initial CRAN submission of mSigSpectra 0.1.1.

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
https://github.com/steverozen/mSigSpectra/actions/runs/25525958885),
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

The "New submission" line is expected — this is the first CRAN
submission of mSigSpectra. The "No suitable spell-checker program
found" line is emitted because the GH Actions runners do not have
aspell installed; it is not a defect of the package.

## Notes seen on local `R CMD check --as-cran`

The local check (with aspell installed) additionally reports:

```
Possibly misspelled words in DESCRIPTION:
  DBS (29:59, 30:42)
  Rozen (34:32)
  SBS (29:54, 30:26)
  VCF (3:48)
  VCFs (27:40)
  al (34:41)
  et (34:38)
  indel (29:64)
  transcriptional (28:59)
```

All of these are intentional: `DBS`, `SBS`, `VCF`, `VCFs`, `indel`,
and `transcriptional` are standard domain terms in mutational-
signature analysis; `Rozen` is the maintainer's surname; `et al` is
part of a Nik-Zainal et al. citation in the Description field.

## Reverse dependencies

There are no reverse dependencies on CRAN (this is a new package).
