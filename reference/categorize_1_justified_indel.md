# Given a indel and its sequence context, categorize it

This function is primarily for internal use, but we export it to
document the underlying logic.

## Usage

``` r
categorize_1_justified_indel(context, ins_or_del, ins_or_del_seq, pos)
```

## Arguments

- context:

  The sequence surrounding the indel PRIOR to the insertion or deletion.

- ins_or_del:

  A singgle character, with "i" denoting an insertion and "d" denotine
  an deletion.

- ins_or_del_seq:

  The the sequence that was inserted or deleted.

- pos:

  For deletions, the 1-based position of the start of the deleted
  sequence; for insertions, the position immediately to the right of
  where the insertion occurs.

## Value

A string that is the canonical representation of the given deletion
type. Return `NA` and raise a warning if there is an un-normalized
representation of the deletion of a repeat unit. See `FindDelMH` for
details. (This seems to be very rare.)

## Details

See
<https://github.com/steverozen/ICAMS/blob/v3.0.9-branch/data-raw/PCAWG7_indel_classification_2021_09_03.xlsx>
for additional information on deletion mutation classification.

This function first handles deletions in homopolymers, then handles
deletions in simple repeats with longer repeat units (e.g. `CACACACA`),
and if the deletion is not in a simple repeat, looks for microhomology.

## Examples

``` r
simplify = function(ll) unlist(ll[c("COSMIC_83", "Koh_89", "Koh_476")])
categorize_1_justified_indel("GGAAAGG", "d", ins_or_del_seq = "A", pos = 3) # "DEL:T:1:2"
#> $ins_or_del
#> [1] "d"
#> 
#> $pre
#> [1] "C"
#> 
#> $ins_or_del_seq
#> [1] "T"
#> 
#> $post
#> [1] "C"
#> 
#> $L
#> [1] 1
#> 
#> $U_seq
#> [1] "T"
#> 
#> $U
#> [1] 1
#> 
#> $U_seq_count_in_indel_seq
#> [1] 1
#> 
#> $indel_str_count_in_ref
#> [1] 3
#> 
#> $R
#> [1] 3
#> 
#> $R_outside_ins_or_del_seq
#> [1] 2
#> 
#> $mh
#> [1] 0
#> 
#> $short_visual
#> [1] "<A>[AA]"
#> 
#> $long_visual
#> [1] "GG <A>[AA] GG"
#> 
#> $unit
#> [1] "A"
#> 
#> $unit_length
#> [1] 1
#> 
#> $internal_rep
#> [1] ""
#> 
#> $internal_reps
#> [1] 0
#> 
#> $spacer
#> [1] ""
#> 
#> $spacer_length
#> [1] 0
#> 
#> $prime3_rep
#> [1] "AA"
#> 
#> $prime3_reps
#> [1] 2
#> 
#> $original_reps
#> [1] 3
#> 
#> $COSMIC_83
#> [1] "DEL:T:1:2"
#> 
#> $Koh_89
#> [1] "C[Del(T):R(1,4)]C"
#> 
#> $Koh_476
#> [1] "C[Del(T):R3]C"
#> 
categorize_1_justified_indel("GGAAAGG", "d", ins_or_del_seq = "A", pos = 4) # "DEL:T:1:2"
#> $ins_or_del
#> [1] "d"
#> 
#> $pre
#> [1] "C"
#> 
#> $ins_or_del_seq
#> [1] "T"
#> 
#> $post
#> [1] "T"
#> 
#> $L
#> [1] 1
#> 
#> $U_seq
#> [1] "T"
#> 
#> $U
#> [1] 1
#> 
#> $U_seq_count_in_indel_seq
#> [1] 1
#> 
#> $indel_str_count_in_ref
#> [1] 2
#> 
#> $R
#> [1] 2
#> 
#> $R_outside_ins_or_del_seq
#> [1] 2
#> 
#> $mh
#> [1] 0
#> 
#> $short_visual
#> [1] "<A>[A]"
#> 
#> $long_visual
#> [1] "GGA <A>[A] GG"
#> 
#> $unit
#> [1] "A"
#> 
#> $unit_length
#> [1] 1
#> 
#> $internal_rep
#> [1] ""
#> 
#> $internal_reps
#> [1] 0
#> 
#> $spacer
#> [1] ""
#> 
#> $spacer_length
#> [1] 0
#> 
#> $prime3_rep
#> [1] "A"
#> 
#> $prime3_reps
#> [1] 1
#> 
#> $original_reps
#> [1] 2
#> 
#> $COSMIC_83
#> [1] "DEL:T:1:1"
#> 
#> $Koh_89
#> [1] "C[Del(T):R(1,4)]T"
#> 
#> $Koh_476
#> [1] "C[Del(T):R2]T"
#> 
simplify(
  categorize_1_justified_indel("TTATT", "d", ins_or_del_seq = "A", pos = 3))
#>           COSMIC_83              Koh_89             Koh_476 
#>         "DEL:T:1:0" "A[Del(T):R(1,4)]A"     "A[Del(T):R1]A" 
simplify(
  categorize_1_justified_indel("TTATATAT", "d", ins_or_del_seq = "TATA", pos = 2))
#> Called from: categorize_1_justified_indel("TTATATAT", "d", ins_or_del_seq = "TATA", 
#>     pos = 2)
#> debug: retlist$COSMIC_83 = gen_COSMIC_83_string(retlist)
#> debug: retlist$Koh_89 = gen_Koh_89_string(retlist)
#> debug: if (retlist$Koh_89 == "Del(2,8):U(1,2):R(2,4)") {
#>     stopifnot(R <= 8)
#> }
#> debug: retlist$Koh_476 = gen_Koh_476_string(retlist)
#> debug: return(retlist)
#>         COSMIC_83            Koh_89           Koh_476 
#>      "DEL:MH:4:3" "Del(4,5):M(3,4)"         "Del4:M3" 
simplify(
  categorize_1_justified_indel("TTATATAT", "d", ins_or_del_seq = "TATAT", pos = 2))
#> Called from: categorize_1_justified_indel("TTATATAT", "d", ins_or_del_seq = "TATAT", 
#>     pos = 2)
#> debug: retlist$COSMIC_83 = gen_COSMIC_83_string(retlist)
#> debug: retlist$Koh_89 = gen_Koh_89_string(retlist)
#> debug: if (retlist$Koh_89 == "Del(2,8):U(1,2):R(2,4)") {
#>     stopifnot(R <= 8)
#> }
#> debug: retlist$Koh_476 = gen_Koh_476_string(retlist)
#> debug: return(retlist)
#>          COSMIC_83             Koh_89            Koh_476 
#> "DEL:repeats:5+:0"       "Del(5,):R1"    "Del5:U(2,):R1" 
```
