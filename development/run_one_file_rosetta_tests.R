# Run build_one_file_rosetta() on each BestMatch476_1 sample in
# ../unified_indels_data/connection_table.tsv and write one .xlsx per sample
# to development/one_file_rosetta_output/.
devtools::load_all(quiet = TRUE)
conn <- read.delim(
  here::here("..", "unified_indels_data", "connection_table.tsv"),
  stringsAsFactors = FALSE
)
ids <- sort(unique(conn$BestMatch476_1))
out_dir <- here::here("development", "one_file_rosetta_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
base <- "/home/steve/MEGA/important_mut_sig_data"
log <- character(0)
for (id in ids) {
  sub <- if (startsWith(id, "SP")) "pcawg_indel_vcfs" else "fmh-unfiltered_vcfs"
  vcf <- file.path(base, sub, paste0(id, ".annotated.indel.vcf.gz"))
  if (!file.exists(vcf)) {
    log <- c(log, paste(id, "MISSING", vcf))
    next
  }
  out <- file.path(out_dir, paste0(id, "_rosetta.xlsx"))
  res <- tryCatch(
    {
      d <- build_one_file_rosetta(vcf, out)
      paste(id, "OK rows =", nrow(d))
    },
    warning = function(w) paste(id, "WARNING", conditionMessage(w)),
    error = function(e) paste(id, "ERROR", conditionMessage(e))
  )
  log <- c(log, res)
}
writeLines(log, file.path(out_dir, "run_log.txt"))
cat(log, sep = "\n")
