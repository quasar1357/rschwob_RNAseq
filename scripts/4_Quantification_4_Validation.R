
quant_path <- "../4_Quantification_3_ExprLevels"

samples <- c("1_1", "1_2", "1_5", "P1", "P2", "P3")

#Sum tmp values from abundance tables for validation (should sum up to 1'000'000)

get_tot_tpm <- function(sample){
  sample_path <- paste0(quant_path, "/", sample, "/abundance.tsv")
  abundance <- read.table(file = sample_path, header = TRUE)
  tot_tpm <- sum(abundance$tpm)
  return(tot_tpm)
}


output <- ""
for (sample in samples){
  output <- paste0(output, sample, ": ", get_tot_tpm(sample), "\n")
}

writeLines(output, "4_Quantification_4_Validation_RESULTS.txt")
