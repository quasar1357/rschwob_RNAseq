# load packages
library(tidyverse)
library(rtracklayer)


############ Load the merged assembly gtf and put in data frame

merged_assembly_gtf <- rtracklayer::import("../results/3_Assembly_2_MergeAssemblies/merged_assembly.gtf")
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
wt_table <- read.csv(file = "../results/5_DiffExpr_3_DiffExpr_TransLevel/sleuth_wt_transcript.csv", header = TRUE)
# lrt_table <- read.csv(file = "../results/5_DiffExpr_3_DiffExpr_TransLevel/sleuth_lrt_transcript.csv", header = TRUE)

log2_and_qval <- function(id, test_table = wt_table){
  entry <- test_table[test_table$target_id == id,]
  log2 <- entry$b
  qval <- entry$qval
  output <- c(log2, qval)
  return(output)
}


### Find the protein coding potential
prot_coding_table <- as.data.frame(read.table("../results/6_IntAn_3_ProtCodPot/cpat_out", header = TRUE, sep="\t"))
prot_coding_pot <- function(id){
  return(prot_coding_table[id,]$coding_prob)
}


### Check if the transcript has a TSS, polyA tail and if it is intergenic
TSS_table <- as.data.frame(read.table("../results/6_IntAn_2_Find_TSS_PolyA_intergenic/TSS_intersect.tsv", header = FALSE, sep="\t"))
TSS_transcripts <- TSS_table[,4]
  
polyA_table <- as.data.frame(read.table("../results/6_IntAn_2_Find_TSS_PolyA_intergenic/polyA_intersect.tsv", header = FALSE, sep="\t"))
polyA_transcripts <- polyA_table[,4]
  
intergenic_table <- as.data.frame(read.table("../results/6_IntAn_2_Find_TSS_PolyA_intergenic/total_intersect.tsv", header = FALSE, sep="\t"))
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
results <- results %>% select(transcript_id, gene_id, gene_name, width)
type <- case_when(is.na(results$gene_name) ~"novel", !is.na(results$gene_name) ~"annotated")
results <- data.frame(results[,1:3], type, results[,4:ncol(results)])

# Check for uniquness of each entry of transcript ID
length(results$transcript_id) == length(unique(results$transcript_id))

# Add columns to be filled later
results[,c("Num_Exons", "prot_coding_pot", "TSS", "PolyA", "Intergenic", "log2_fold_change", "q_val")] <- NA
num_transcripts <- length(results$transcript_id)

# Look up and add the values for all transcripts
# NOTE: this can take a long time to compute and can be done step by step for the different attributes
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


############ Write, filter (and read) the results

### Write the results into a file (or read from a saved file, so it doesn't have to compute everything again)
write.csv(results, file = "../results/7_Summary/results_all.csv", row.names = FALSE)
# results <- read.csv(file = "../results/7_Summary/results_all.csv", header = TRUE)

### Filter the results for the most interesting entries, rank them and write them into a file (or read from a saved file)
width_threshold <- 200
log2_fold_threshold <- 2
q_val_threshold <- 0.025 # Because we are testing 2-sided! This mean the significance level is 0.05 
prot_coding_pot_threshold <- 0.364

results_filtered <- results %>% filter(Width >= width_threshold, prot_coding_pot <= prot_coding_pot_threshold, Num_Exons > 1, log2_fold_change < -log2_fold_threshold | log2_fold_change > log2_fold_threshold, q_val <= q_val_threshold, TSS == TRUE, PolyA == TRUE)

# ADD SCORE
norm_abs <- function(vec){
  vec_abs <- abs(vec)
  vec_norm <- (vec_abs - min(vec_abs)) / (max(vec_abs) - min(vec_abs))
}

fold_norm <- norm_abs(results_filtered$log2_fold_change)
q_log <- -log10(results_filtered$q_val)
q_norm <- norm_abs(q_log)
fold_and_q <- (fold_norm + q_norm) / 2
fold_and_q_norm <- norm_abs(fold_and_q)

prot_cod_norm <- norm_abs(results_filtered$prot_coding_pot)

score <- fold_and_q_norm - prot_cod_norm
results_filtered$score <- score

# SORT AND ADD RANK
results_filtered <- results_filtered %>% arrange(desc(score))
results_filtered$rank <- 1:nrow(results_filtered)
results_filtered <- results_filtered  %>% arrange(desc(type))

# For novel and annotated together
# results_filtered <- read.csv(file = "../results/7_Summary/results_filtered.csv", header = TRUE)
write.csv(results_filtered, file = "../results/7_Summary/results_filtered.csv", row.names = FALSE)

# For novel and annotated separately
# results_filtered_novel <- results_filtered %>% filter(is.na(gene_name))
# write.csv(results_filtered_novel, file = "../results/7_Summary/results_filtered_novel.csv", row.names = FALSE)
# results_filtered_annotated <- results_filtered %>% filter(is.na(gene_name) == FALSE)
# write.csv(results_filtered_annotated, file = "../results/7_Summary/results_filtered_annotated.csv", row.names = FALSE)


### Filter and save GTF for those transcripts (for further use...)
# merged_assembly_gtf_filtered <- merged_assembly_gtf[which(merged_assembly_gtf$transcript_id %in% results_filtered$transcript_id),]
# merged_assembly_filtered_df <- as.data.frame(merged_assembly_filtered)
# rtracklayer::export(merged_assembly_gtf_filtered, "../results/7_Summary/merged_assembly_filtered.gtf")

