# Map a character ref_genome argument to its canonical BSgenome package name

Map a character ref_genome argument to its canonical BSgenome package
name

## Usage

``` r
infer_ref_genome_name(ref_genome)
```

## Value

A single character string giving the canonical 'BSgenome' package name,
e.g. `"BSgenome.Hsapiens.UCSC.hg38"`. Errors if `ref_genome` is not
recognized.
