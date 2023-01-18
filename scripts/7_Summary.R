# load packages
library(tidyverse)
library(rtracklayer)


############ Load the merged assembly gtf and put in data frame

merged_assembly_gtf <- rtracklayer::import("../3_Assembly_2_MergeAssemblies/merged_assembly.gtf")
merged_assembly <- as.data.frame(merged_assembly_gtf)


############ Create functions to find the interesting values

### Find the number of exons from the merged assembly gtf
merged_assembly_exons <- merged_assembly %>% select(transcript_id, exon_number)
merged_assembly_exons_sorted <- merged_assembly_exons[order(merged_assembly_exons$exon_number, decreasing = TRUE),]

num_exons <- function(id){
  return(merged_assembly_exons_sorted$exon_number[match(id, merged_assembly_exons_sorted$transcript_id)])
}

# num_exons_slow <- function(id){
#   exons <- merged_assembly %>% filter(type == "exon", transcript_id == id)
#   return(nrow(exons))
# }


### Find the log2-fold change and the q-value from the Differential Expression tables (on transcript level)
wt_table <- read.csv(file = "../5_DiffExpr_3_DifExpr_TransLevel/sleuth_wt_transcript.csv", header = TRUE)
# lrt_table <- read.csv(file = "../5_DiffExpr_3_DifExpr_TransLevel/sleuth_lrt_transcript.csv", header = TRUE)

log2_and_qval <- function(id, test_table = wt_table){
  entry <- test_table[test_table$target_id == id,]
  log2 <- entry$b
  qval <- entry$qval
  output <- c(log2, qval)
  return(output)
}


### Find the protein coding potential
prot_coding_table <- as.data.frame(read.table("../6_IntAn_3_ProtCodPot/cpat_out", header = TRUE, sep="\t"))
prot_coding_pot <- function(id){
  return(prot_coding_table[id,]$coding_prob)
}


### Check if the transcript has a TSS, polyA tail or if it is intergenic
TSS_table <- as.data.frame(read.table("../6_IntAn_2_Find_TSS_PolyA_intergenic/TSS_intersect.tsv", header = FALSE, sep="\t"))
TSS_transcripts <- TSS_table[,4]
  
polyA_table <- as.data.frame(read.table("../6_IntAn_2_Find_TSS_PolyA_intergenic/polyA_intersect.tsv", header = FALSE, sep="\t"))
polyA_transcripts <- polyA_table[,4]
  
intergenic_table <- as.data.frame(read.table("../6_IntAn_2_Find_TSS_PolyA_intergenic/total_intersect.tsv", header = FALSE, sep="\t"))
intergenic_transcripts <- intergenic_table[,4]

has_TSS <- function(id){
  return(id %in% TSS_transcripts)
}

has_polyA <- function(id){
  return(id %in% polyA_transcripts)
}

is_intergenic <- function(id){
  return(id %in% intergenic_transcripts)
}


############ Run and use the functions to fill in the infos into the results data frame

# For the results, we filter the merged assembly for all transcripts and take only the columns we're interested in
results <- merged_assembly %>% filter(type == "transcript")
results <- results %>% select(transcript_id, gene_id, gene_name)

# Check for uniquness of each entry of transcript ID
length(results$transcript_id) == length(unique(results$transcript_id))

# Add columns to be filled later
results[,c("Num_Exons", "log2_fold_change", "q_val", "prot_coding_pot", "TSS", "PolyA", "Intergenic")] <- NA
num_transcripts <- length(results$transcript_id)

# Look up and add the values for all transcripts
for(idx in seq(1, num_transcripts)){
  id <- results$transcript_id[idx]
  
  results$Num_Exons[idx] <- num_exons(id)

  test_results <- log2_and_qval(id)
  results$log2_fold_change[idx] <- test_results[1]
  results$q_val[idx] <- test_results[2]

  results$prot_coding_pot[idx] <- prot_coding_pot(id)
  
  # results$TSS[idx] <- has_TSS(id)
  # results$PolyA[idx] <- has_polyA(id)
  # results$Intergenic[idx] <- is_intergenic(id)
}

# For TSS, PolyA and intergenic we only need to look up all transcripts that are in those tables
# All the other transcripts are assumed to not have/be TSS, PolyA or intergenic
results$TSS <- FALSE
for(t in TSS_transcripts){
  results$TSS[match(t, results$transcript_id)] <- TRUE
}
sum(results$TSS) == length(TSS_transcripts)

results$PolyA <- FALSE
for(t in polyA_transcripts){
  results$PolyA[match(t, results$transcript_id)] <- TRUE
}
sum(results$PolyA) == length(polyA_transcripts)

results$Intergenic <- FALSE
for(t in intergenic_transcripts){
  results$Intergenic[match(t, results$transcript_id)] <- TRUE
}
sum(results$Intergenic) == length(intergenic_transcripts)


############ Write, filter and read the results

### Write the results into a file (or read from a saved file, so it doesn't have to compute everything again)
write.csv(results, file = "../7_Summary/results_all.csv", row.names = FALSE)
# results <- read.csv(file = "../7_Summary/results_all.csv", header = TRUE)

### Filter the results for the most interesting entries and write them into a file (or read from a saved file)
log2_fold_threshold <- 2
q_val_threshold <- 0.05
prot_coding_pot_threshold <- 0.364

# For novel and annotated together
results_filtered <- results %>% filter(prot_coding_pot <= prot_coding_pot_threshold, Num_Exons > 1, log2_fold_change < -log2_fold_threshold | log2_fold_change > log2_fold_threshold, q_val <= q_val_threshold, TSS == TRUE, PolyA == TRUE)
write.csv(results_filtered, file = "../7_Summary/results_filtered.csv", row.names = FALSE)
# results_filtered <- read.csv(file = "../7_Summary/results_filtered.csv", header = TRUE)
coding <- sum(results$prot_coding_pot > prot_coding_pot_threshold)
non_coding <- sum(results$prot_coding_pot <= prot_coding_pot_threshold)
total <- coding + non_coding
coding_perc <- 100 * coding / total
non_coding_perc <- 100 * non_coding / total

# For novel and annotated separately
results_filtered_novel <- results %>% filter(is.na(gene_name), prot_coding_pot < prot_coding_pot_threshold, Num_Exons > 1, log2_fold_change < -log2_fold_threshold | log2_fold_change > log2_fold_threshold, q_val <= q_val_threshold, TSS == TRUE, PolyA == TRUE)
write.csv(results_filtered_novel, file = "../7_Summary/results_filtered_novel.csv", row.names = FALSE)
results_filtered_annotated <- results %>% filter(is.na(gene_name) == FALSE, prot_coding_pot <= prot_coding_pot_threshold, Num_Exons > 1, log2_fold_change < -log2_fold_threshold | log2_fold_change > log2_fold_threshold, q_val <= q_val_threshold, TSS == TRUE, PolyA == TRUE)
write.csv(results_filtered_annotated, file = "../7_Summary/results_filtered_annotated.csv", row.names = FALSE)


### Filter and save GTF for those transcripts (for further use...)
merged_assembly_gtf_filtered <- merged_assembly_gtf[which(merged_assembly_gtf$transcript_id %in% results_filtered$transcript_id),]
# merged_assembly_filtered_df <- as.data.frame(merged_assembly_filtered)
rtracklayer::export(merged_assembly_gtf_filtered, "../7_Summary/merged_assembly_filtered.gtf")
