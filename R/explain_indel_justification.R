explain_indel_justification = function(
  long_str,
  short_str,
  pos_before,
  pos_after,
  is_ins
) {
  len_diff = nchar(long_str) - nchar(short_str)

  message(
    "\n\nJustification explanation: ",
    ifelse(is_ins, "insertion", "deletion"),
    " of ",
    stringi::stri_sub(
      long_str,
      pos_before,
      pos_before + len_diff - 1
    ),
    " ========="
  )
  message("Prior to justifying: pos = ", pos_before)
  message("After justifying:   pos = ", pos_after)
  message(
    "After justifying the ",
    ifelse(is_ins, "inserted", "deleted"),
    " sequence is ",
    stringi::stri_sub(
      long_str,
      pos_after,
      pos_after + len_diff - 1
    )
  )
}
