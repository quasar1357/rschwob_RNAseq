# load packages
library(tidyverse)
library(rtracklayer)

results <- read.csv(file = "../7_Summary/results_all.csv", header = TRUE)

width_threshold <- 200
log2_fold_threshold <- 2
q_val_threshold <- 0.025 # Because we are testing 2-sided! Like this, the significance level is 0.05 
prot_coding_pot_threshold <- 0.364


########## Get the proportion of all important statistics (separately for each value)

total <- nrow(results)
total_quant <- total - sum(is.na(results$log2_fold_change))

# Manually check that there are only 2 values (TRUE and FALSE) for TSS, polyA and intergenic --> can only count one and divide by total
# print(length(unique(results[,10])))

### Get proportion novel vs. annotated
novel <- sum(is.na(results$gene_name))
prop_novel <- novel / total

### Get proportion of long (>= 200nt) transcripts
long <- sum(results$Width >= width_threshold)
prop_long <- long / total

### Get proportion multi-exon vs. single-exon
multi_exon <- sum(results$Num_Exons > 1)
prop_multi_exon <- multi_exon / total
# single_exon <- sum(results$Num_Exons == 1)

### Get proportion noncoding vs. coding
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

### Get proportion outside the log2_fold_change threshold (if we choose 2, it would mean: -2 > log_2_fold > 2)
# NOTE: here, we want to exclude the NAs, so the transcripts that could not have been evaluated in the analysis
high_fold <- sum(results$log2_fold_change < -log2_fold_threshold | results$log2_fold_change > log2_fold_threshold, na.rm = TRUE)
prop_high_fold <- high_fold / total
prop_high_fold_per_quant <- high_fold / total_quant

### Get proportion significant vs. non-significant (here with q_val <= 0.05)
# NOTE: here, we want to exclude the NAs, so the transcripts that could not have been evaluated in the analysis
significant <- sum(results$q_val <= q_val_threshold, na.rm = TRUE)
prop_significant <- significant / total
prop_significant_per_quant <- significant / total_quant


### SUMMARIZE AND SAVE ALL THOSE STATISTICS
results_stats <- data.frame(prop_novel, prop_long, prop_multi_exon, prop_non_coding, prop_TSS, prop_PolyA, prop_Intergenic, prop_high_fold, prop_significant, prop_high_fold_per_quant, prop_significant_per_quant)
write.csv(results_stats, "../7_Summary/results_stats.csv", row.names = FALSE)

results_stats_perc <- results_stats * 100


########## Create functions to analyze different combinations of those attributes

### Functions used to quantify
num <- function(df){
  return(nrow(df))
}

prop <- function(df){
  return(num(df) / total)
}

perc <- function(df){
  return(100 * prop(df))
}

# Use these functions to find the proportion/percent among entries that have a valid quantification result (i.e. no NA in the columns for log2_fold_change and q_val)
prop_quant <- function(df){
  return(num(df) / total_quant)
}

perc_quant <- function(df){
  return(100 * prop_quant(df))
}

# This function serves to check which/how many of the chosen entries have a valid quantification result (i.e. no NA in the columns for log2_fold_change and q_val)
has_quant <- function(df){
  filtered <- df %>% filter(is.na(log2_fold_change) == FALSE, is.na(q_val) == FALSE)
return(filtered)
}


### Functions used to filter (potentially consecutively...)
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

# This is a wrapper function to find solid lncRNA candidates that have all the necessary prerequisites
solid_candidates <- function(df){
  return(long(multi_exon(non_coding(tss(polya(df))))))
}

high_fold <- function(df){
  filtered <- df %>% filter(log2_fold_change < -log2_fold_threshold | log2_fold_change > log2_fold_threshold, na.rm = TRUE)
  return(filtered)
}

significant <- function(df){
  filtered <- df %>% filter(q_val <= q_val_threshold, na.rm = TRUE)
  return(filtered)
}

# Again this is a wrapper function to find transcripts with both a high log2_fold_change AND a high significance
sign_high <- function(df){
  return(significant(high_fold(df)))
}


########## Find numbers for combinations of interest

# NOTE: if we use high_fold or significant with prop or perc, we want to use the corrected total that excludes entries that have NA in those columns
total <- nrow(results)
total_quant <- total - sum(is.na(results$log2_fold_change))

### Total solid candidate lncRNAs
num(solid_candidates(results))
# 6496
perc(solid_candidates(results))
# 2.933831

### How many were even correctly quantified in total?
num(has_quant(results))
# 78930
perc(has_quant(results))
# 35.64767

# That makes how many correctly quantified solid candidates?
num(has_quant(solid_candidates(results)))
# 3010
perc(has_quant(solid_candidates(results)))
# 1.359426
# And how many solid candidates AMONG THE CORRECTLY QUANTIFIED?
perc_quant(has_quant(solid_candidates(results)))
# 3.813506

### How many transcripts were significantly expressed with a high log2_fold_change?
num(sign_high(results))
# 1187
perc(sign_high(results))
# 0.5360925
# And how many significant and high_fold AMONG THE CORRECTLY QUANTIFIED?
perc_quant(sign_high(results))
# 1.503864

### RESULT: How many solid candidates are significant and high_fold? --> Those are the "hot" candidates!
num(sign_high(solid_candidates(results)))
# 43
perc(sign_high(solid_candidates(results)))
# 0.01942037
perc_quant(sign_high(solid_candidates(results)))
# 0.05447865

### Now check for novel transcripts
num(novel(results))
# 13106
perc(novel(results))
# 5.919148
num(novel(solid_candidates(results)))
# 433
perc(novel(solid_candidates(results)))
# 0.1955586
num(novel(sign_high(solid_candidates(results))))
# 14
perc(novel(sign_high(solid_candidates(results))))
# 0.006322911
perc_quant(novel(sign_high(solid_candidates(results))))
# 0.01773724

