# load packages
library(tidyverse)
library(rtracklayer)

results <- read.csv(file = "../7_Summary/results_all.csv", header = TRUE)

width_threshold <- 200
log2_fold_threshold <- 2
q_val_threshold <- 0.025 # Because we are testing 2-sided! This mean the significance level is 0.05 
prot_coding_pot_threshold <- 0.364

total <- nrow(results)

# Manually check that there are only 2 values (TRUE and FALSE) for TSS, polyA and intergenic --> can only count one and divide by total
# print(length(unique(results[,10])))

### Get proportion novel vs. annotated
novel <- sum(is.na(results$gene_name))
prop_novel <- novel / total

### Get proportion of long (>= 200nt) transcripts
long <- sum(results$Width >= width_threshold)
prop_long <- long / total

### Get proportion single-exon vs. multi-exon
multi_exon <- sum(results$Num_Exons > 1)
prop_multi_exon <- multi_exon / total
# single_exon <- sum(results$Num_Exons == 1)

### Get proportion outside the log2_fold_change threshold (if we choose 2, it would mean: -2 > log_2_fold > 2)
# NOTE: here, we want to exclude the NAs, so the transcripts that could not have been evaluated in the analysis
high_fold <- sum(results$log2_fold_change < -log2_fold_threshold | results$log2_fold_change > log2_fold_threshold, na.rm = TRUE)
low_fold <- sum(results$log2_fold_change >= -log2_fold_threshold & results$log2_fold_change <= log2_fold_threshold, na.rm = TRUE)
total_evaluated_log2_fold <- high_fold + low_fold
prop_high_fold <- high_fold / total_evaluated_log2_fold

### Get proportion significant vs. non-significant (here with q_val <= 0.05)
# NOTE: here, we want to exclude the NAs, so the transcripts that could not have been evaluated in the analysis
significant <- sum(results$q_val <= q_val_threshold, na.rm = TRUE)
non_significant <- sum(results$q_val > q_val_threshold, na.rm = TRUE)
total_evaluated_qval <- significant + non_significant
prop_significant <- significant / total_evaluated_qval

### Get proportion coding vs. noncoding
non_coding <- sum(results$prot_coding_pot <= prot_coding_pot_threshold)
prop_non_coding <- non_coding / total

# Get proportion TSS
TSS <- sum(results$TSS)
prop_TSS <- TSS / total

# Get proportion polyA
PolyA <- sum(results$PolyA)
prop_PolyA <- PolyA / total

# Get proportion "intergenic"
Intergenic <- sum(results$Intergenic)
prop_Intergenic <- Intergenic / total

### SUMMARIZE AND SAVE ALL THOSE STATISTICS
results_stats <- data.frame(prop_novel, prop_long, prop_multi_exon, prop_high_fold, prop_significant, prop_non_coding, prop_TSS, prop_PolyA, prop_Intergenic)
write.csv(results_stats, "../7_Summary/result_stats.csv", row.names = FALSE)

results_stats_perc <- results_stats * 100


########## Create functions to analyze different combinations of those attributes

# Functions used to quantify
num <- function(df){
  return(nrow(df))
}

prop <- function(df){
  return(num(df) / total)
}

perc <- function(df){
  return(100 * prop(df))
}

# Functions used to filter (potentially consecutively...)
novel <- function(df){
  filtered <- df %>% filter(is.na(gene_name))
  return(filtered)
}

long <- function(df){
  filtered <- df %>% filter(Width >= width_threshold)
  return(filtered)
}

multi_exon <- function(df){
  filtered <- df %>% filter(Num_Exons > 1)
  return(filtered)
}

high_fold <- function(df){
  filtered <- df %>% filter(log2_fold_change < -log2_fold_threshold | log2_fold_change > log2_fold_threshold, na.rm = TRUE)
  return(filtered)
}

significant <- function(df){
  filtered <- df %>% filter(q_val <= q_val_threshold, na.rm = TRUE)
  return(filtered)
}

non_coding <- function(df){
  filtered <- df %>% filter(prot_coding_pot <= prot_coding_pot_threshold)
  return(filtered)
}

tss <- function(df){
  filtered <- df %>% filter(TSS)
  return(filtered)
}

polya <- function(df){
  filtered <- df %>% filter(PolyA)
  return(filtered)
}

intergenic <- function(df){
  filtered <- df %>% filter(Intergenic)
  return(filtered)
}

# NOTE: if you use high_fold or significant with prop or perc, you want to use the corrected total that excludes entries that have NA in those columns!
total <- nrow(results)
total <- total - sum(is.na(results$log2_fold_change))

test <- long(non_coding(tss(polya(results))))


