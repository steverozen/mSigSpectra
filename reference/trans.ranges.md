# Transcript ranges for transcriptional strand annotation

Precomputed transcript ranges (one row per gene) used by
[`add_transcript_strand()`](https://steverozen.github.io/mSigSpectra/reference/add_transcript_strand.md)
to determine the coding strand for each variant. Sources: GENCODE v30
(human) and vM21 (mouse). Only genes with CCDS IDs are retained.

## Usage

``` r
trans.ranges.GRCh37

trans.ranges.GRCh38

trans.ranges.GRCm38
```

## Format

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `chrom`, `start`, `end`, `strand`, `Ensembl.gene.ID`,
`gene.symbol`. One-based coordinates.

An object of class `data.table` (inherits from `data.frame`) with 19083
rows and 6 columns.

An object of class `data.table` (inherits from `data.frame`) with 19096
rows and 6 columns.

An object of class `data.table` (inherits from `data.frame`) with 20325
rows and 6 columns.

## Source

<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gff3.gz>

<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gff3.gz>

<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz>
