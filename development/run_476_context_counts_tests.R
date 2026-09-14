# Test build_476_context_counts() on DRUP01030028T and SP48008 using each
# sample's 476-types Del2:U2:R2 and Del2:M1. Writes to
# development/one_file_rosetta_output/<id>_context_counts.xlsx.
devtools::load_all(quiet = TRUE)
source(here::here("development", "build_476_context_counts.R"))
base <- "/home/steve/MEGA/important_mut_sig_data"
out_dir <- here::here("development", "one_file_rosetta_output")
vcfs <- c(
  DRUP01030028T = file.path(
    base, "fmh-unfiltered_vcfs", "DRUP01030028T.annotated.indel.vcf.gz"
  ),
  SP48008 = file.path(
    base, "pcawg_indel_vcfs", "SP48008.annotated.indel.vcf.gz"
  )
)
types <- list(
  DRUP01030028T = c("Del2:U2:R2", "Del2:M1"),
  SP48008 = c("Del2:U2:R2", "Del2:M1")
)
for (id in names(vcfs)) {
  out <- file.path(out_dir, paste0(id, "_context_counts.xlsx"))
  d <- build_476_context_counts(
    vcfs[[id]], types[[id]], out, left_num = 2, right_num = 5
  )
  cat("\n==", id, "rows =", nrow(d), "\n")
  print(d[, utils::head(.SD, 4), by = Koh_476][, .(Koh_476, n_indels, context, n)])
}
