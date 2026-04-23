# Add sequence context and transcript information to an in-memory ID (insertion/deletion) VCF, and confirm that they match the given reference genome

Add sequence context and transcript information to an in-memory ID
(insertion/deletion) VCF, and confirm that they match the given
reference genome

## Usage

``` r
justify_id_vcf(
  ID.vcf,
  ref.genome,
  name.of.VCF = NULL,
  suppress.discarded.variants.warnings = TRUE,
  explain_indels = 1,
  context_width_multiplier = 20L
)
```

## Arguments

- ID.vcf:

  An in-memory ID (insertion/deletion) VCF as a `data.frame`. This
  function expects that there is a "context base" to the left, for
  example REF = ACG, ALT = A (deletion of CG) or REF = A, ALT = ACC
  (insertion of CC).

- ref.genome:

  Can be a string or a BSgenome. If a string, it should a well-known
  name for reference genome in BSgenome

- name.of.VCF:

  Name of the VCF file.

- suppress.discarded.variants.warnings:

  If TRUE, do warn when variants that cannot be processed are discarded.

- explain_indels:

  If 0, do not explain, if 1, explain indels that were justifed (via
  messages), if 2, geenrate messages for all indels.

- context_width_multiplier:

  Used to guess how much sequence on each side of an indel is needed to
  categorize it.

## Value

A list of elements:

- `annotated.vcf`: The original VCF data frame in which `POS`, `REF` and
  `ALT` with new columns added to the input data frame:

  - `seq.context`: The sequence embedding the variant.

  - `seq.context.width`: The width of `seq.context` to the left of the
    indel

  - `pos_shift`

- `discarded.variants`: **Non-NULL only if** there are variants that
  were excluded from the analysis. See the added extra column
  `discarded.reason` for more details.

## Examples

``` r
# See tests/testthat/test_indel_classification.R for an end-to-end example
# against a real Strelka ID VCF.
```
