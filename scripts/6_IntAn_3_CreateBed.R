library(tidyverse)
library(rtracklayer)

############ Function to load a gtf and convert it to a data frame

load_df_gtf <- function(gtf_path){
  # Load GTF and convert it to a data frame
  gtf_in <- rtracklayer::import(gtf_path)
  gtf_df <- as.data.frame(gtf_in)
  
  # Optional: Only take important columns
  gtf <- gtf_df # %>% select(seqnames, strand, type, gene_id, transcript_id, exon_number, gene_name, ref_gene_id)
  return(gtf)
}


gtf <- load_df_gtf("../3_Assembly_2_MergeAssemblies/merged_assembly.gtf")

# Filter the gtf data frame for transcripts that were found on a chromosome, select columns for the bed format
gtf_filtered <- gtf %>% filter(type == "transcript") %>% filter(str_detect(seqnames, "chr")) %>% select(seqnames, start, end, transcript_id, score, strand)
gtf_for_bed <- gtf_filtered

names(gtf_for_bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

gtf_for_bed[,c("thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")] <- 0
gtf_for_bed$thickStart <- gtf_for_bed$chromStart
gtf_for_bed$thickEnd <- gtf_for_bed$chromEnd
gtf_for_bed$itemRgb <- "255,0,0"

write.table(gtf_for_bed, file = "../6_IntAn_3_ProtCodPot/merged_assembly.bed", row.names = F, col.names = F, quote = F, sep = "\t")
