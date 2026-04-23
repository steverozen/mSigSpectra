# K-mer abundances for density calculations

Nested list of k-mer counts keyed by
`BSgenome.Hsapiens.1000genomes.hs37d5` / `BSgenome.Hsapiens.UCSC.hg38` /
`BSgenome.Mmusculus.UCSC.mm10`, then by `exome` / `transcript` /
`genome`, then by catalog size as a string (`"78"`, `"96"`, `"136"`,
`"144"`, `"192"`, `"1536"`). Each leaf value is a named integer vector:
k-mer → count.

## Usage

``` r
all.abundance
```

## Format

An object of class `list` of length 3.

## Details

Used to convert counts to density (mutations per megabase of context).
See the planned
[`transform_catalog()`](https://steverozen.github.io/mSigSpectra/reference/transform_catalog.md)
function (not yet ported).

## Examples

``` r
all.abundance$BSgenome.Hsapiens.UCSC.hg38$transcript$`144`[1:5]
#>       AA       AC       AG       AT       CA 
#> 90769160 57156295 85738416 87552737 83479655 
```
