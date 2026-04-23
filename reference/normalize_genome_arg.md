# Normalize a reference-genome argument to a BSgenome object

Accepts either a `BSgenome` object directly, or one of the string
identifiers `"GRCh37"` / `"hg19"` /
`"BSgenome.Hsapiens.1000genomes.hs37d5"`, `"GRCh38"` / `"hg38"` /
`"BSgenome.Hsapiens.UCSC.hg38"`, or `"GRCm38"` / `"mm10"` /
`"BSgenome.Mmusculus.UCSC.mm10"`. The relevant BSgenome package must be
installed; an informative error is raised if not.

## Usage

``` r
normalize_genome_arg(ref_genome)
```

## Arguments

- ref_genome:

  A BSgenome object or a recognized character identifier.

## Value

A BSgenome object.
