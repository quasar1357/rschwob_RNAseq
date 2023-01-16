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

############ Create bed files for the merged assembly gtf

### 1) Preparations & creating a bed file for the total transcripts

gtf <- load_df_gtf("../3_Assembly_2_MergeAssemblies/merged_assembly.gtf")

# Filter the gtf data frame for transcripts that were found on a chromosome, select columns for the bed format
gtf_filtered <- gtf %>% filter(type == "transcript") %>% filter(str_detect(seqnames, "chr")) %>% select(seqnames, start, end, transcript_id, score, strand)
# Write into a file
write.table(gtf_filtered, file = "../6_IntAn_1_CreateBed/merged_assembly_TOTAL.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Get the plus and minus stranded transcripts
is_plus <- which(gtf_filtered$strand == "+")
is_minus <- which(gtf_filtered$strand == "-")


###  2) Create a bed file for the analysis of the TSS.

# We always want the BEGINNING of the transcript: If the transcript is on the "+" strand we need the "start" coordinates, if it is on the "-" strand we need the "end" coordinates.
# We want a 100bp window, we need to add +-50

TSS_plus <- gtf_filtered[is_plus, ]
TSS_plus$end <- TSS_plus$start + 50
TSS_plus$start <- TSS_plus$start - 50

TSS_minus <- gtf_filtered[is_minus, ]
TSS_minus$start <- TSS_minus$end - 50
TSS_minus$end <- TSS_minus$end + 50

# Save the adjusted data for plus and minus strands in one data frame
TSS_bed <- as.data.frame(rbind(TSS_plus, TSS_minus))
# Reorder the data as in the initial gtf file
TSS_bed <- TSS_bed[order(match(TSS_bed$transcript_id, gtf_filtered$transcript_id)), ]
# In case there are any negative entries, they should be changed to 0
TSS_bed$start[which(TSS_bed$start < 0)] <- 0
# Write into a file
write.table(TSS_bed, file = "../6_IntAn_1_CreateBed/merged_assembly_TSS.bed", row.names = F, col.names = F, quote = F, sep = "\t")


### 3) Create a bed file for the analysis of the polyA.

# We always want the END of the transcript: If the transcript is on the "+" strand we need the "end" coordinates, if it is on the "-" strand we need the "start" coordinates.
# We want a 100bp window, we need to add +-50

polyA_plus <- gtf_filtered[is_plus,]
polyA_plus$start <- polyA_plus$end - 50
polyA_plus$end <- polyA_plus$end + 50

polyA_minus <- gtf_filtered[is_minus, ]
polyA_minus$end <- polyA_minus$start + 50
polyA_minus$start <- polyA_minus$start - 50

# Save the adjusted data for plus and minus strands in one data frame
polyA_bed <- as.data.frame(rbind(polyA_plus, polyA_minus))
# Reorder the data as in the initial gtf file
polyA_bed <- polyA_bed[order(match(polyA_bed$transcript_id, gtf_filtered$transcript_id)), ]  
# In case there are any negative entries, they should be changed to 0
polyA_bed$start[which(polyA_bed$start < 0)] <- 0
# Write into a file
write.table(polyA_bed, file = "../6_IntAn_1_CreateBed/merged_assembly_polyA.bed", row.names = F, col.names = F, quote = F, sep = "\t")


############ Adjust/create reference bed files

### 1) Adjust the bed reference for PolyA analysis

# For the polyA analysis we need to adapt the downloaded bed file.
polyA_bed <- as.data.frame(read.table("../references/atlas.clusters.2.0.GRCh38.96.bed", header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))
names(polyA_bed) <- c("seqnames", "start", "end", "transcript_id", "tpm_1", "strand", "sample_percentage", "nr3_protocols", "tpm_2", "cluster_annotation", "information")
polyA_bed_filtered <- polyA_bed %>% filter(str_detect(seqnames, "^[[:digit:]]+") | seqnames == "X" | seqnames == "Y")
polyA_bed_chr <- polyA_bed_filtered 
polyA_bed_chr$seqnames <- sub("^", "chr", polyA_bed_chr$seqnames)
# Write into a file
write.table(polyA_bed_chr, file = "../6_IntAn_1_CreateBed/atlas.clusters.2.0.GRCh38.96_ADAPTED.bed", row.names = F, col.names = F, quote = F, sep = "\t")


### 2) Create reference bed file for the analysis of "intergenic" genes

# Load GTF and convert it to a data frame
ref_gtf <- load_df_gtf("../references/gencode.v21.chr_patch_hapl_scaff.annotation.gtf")

# Filter the gtf data frame for transcripts that were found on a chromosome, select columns for the bed format
ref_gtf_filtered <- ref_gtf %>% filter(type == "transcript") %>% filter(str_detect(seqnames, "chr")) %>% select(seqnames, start, end, transcript_id, score, strand)
# Add an arbitrary score to the bed entries
ref_gtf_filtered$score <- rep(1000, length.out = nrow(ref_gtf_filtered))
# Write into a file
write.table(ref_gtf_filtered, file = "../6_IntAn_1_CreateBed/reference_TOTAL.bed", row.names = F, col.names = F, quote = F, sep = "\t")
