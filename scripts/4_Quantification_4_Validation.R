
quant_path <- "../results/4_Quantification_3_ExprLevels"

samples <- c("1_1", "1_2", "1_5", "P1", "P2", "P3")

out_path <- "../results/4_Quantification_4_Validation"

#Sum tmp values from abundance tables for validation (should sum up to 1'000'000)

get_tot_tpm <- function(sample){
  sample_path <- paste0(quant_path, "/", sample, "/abundance.tsv")
  abundance <- read.table(file = sample_path, header = TRUE)
  tot_tpm <- sum(abundance$tpm)
  return(tot_tpm)
}

test_tpm <- ""
for (sample in samples){
  test_tpm <- paste0(test_tpm, sample, ": ", get_tot_tpm(sample), "\n")
}

writeLines(test_tpm, paste0(out_path, "/4_Quantification_4_Validation_RESULTS.txt"))



# Find the number of entries in tpm that are >0, i.e. expressed transcripts

get_num_expr <- function(sample){
  sample_path <- paste0(quant_path, "/", sample, "/abundance.tsv")
  abundance <- read.table(file = sample_path, header = TRUE)
  num_expr <- sum(abundance$tpm>0)
  return(num_expr)
}

num_expr <- ""
for (sample in samples){
  num_expr <- paste0(num_expr, sample, ": ", get_num_expr(sample), "\n")
}

writeLines(num_expr, paste0(out_path, "/4_Quantification_4_NumExpr_RESULTS.txt"))