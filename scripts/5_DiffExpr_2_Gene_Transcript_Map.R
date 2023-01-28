library(tidyverse)
library(rtracklayer)

# Load experiment information
experiment_info <- read.csv("../results/5_DiffExpr_1_Experiment_Table/experiment_table.csv")

############  Associate transcripts to genes

# Load GTF and convert it to a data frame
gtf <- rtracklayer::import("../results/3_Assembly_2_MergeAssemblies/merged_assembly.gtf")
gtf_df <- as.data.frame(gtf)

# Filter the columns that are needed and arrange to get NAs at the end
gtf_filtered <- gtf_df %>% filter(type == "transcript") %>% select(gene_id, transcript_id, gene_name, ref_gene_id)
gtf_arranged <- gtf_filtered %>% arrange(gene_name)

# Get the row numbers that are NA and that are not NA
is_NA <- is.na(gtf_arranged$gene_name)
is_NA_nr <- which(is_NA == TRUE)
is_not_NA_nr <- which(is_NA == FALSE)

# Get the gene_id and transcript_id if gene_name is NA
is_NA_nr_df <- gtf_arranged[is_NA_nr, c(1, 2)]
# Get the gene_name and transcript_id if gene_name is not NA
is_not_NA_nr_df <- gtf_arranged[is_not_NA_nr, c(3, 2)]
colnames(is_not_NA_nr_df) <- c("gene_id", "transcript_id")

# Put results in one data frame and sort it the same as the gtf_arranged
gene_transcript_map <- as.data.frame(rbind(is_NA_nr_df, is_not_NA_nr_df))
gene_transcript_map <- gene_transcript_map[order(match(gene_transcript_map[,2],gtf_filtered[,2])),]
colnames(gene_transcript_map) <- c("gene_id", "target_id")

write.csv(gene_transcript_map, file = "../results/5_DiffExpr_2_Gene_Transcript_Map/gene_transcript_map.csv", row.names = FALSE)



