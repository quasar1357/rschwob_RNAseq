# load packages
library(tidyverse)
library(rtracklayer)

# Load merged assembly GTF and convert it to a data frame
load_filter_gtf <- function(gtf_path){
  gtf_in <- rtracklayer::import(gtf_path)
  gtf_df <- as.data.frame(gtf_in)
  
  # Create filter of the gtf df, so we can get other numbers easily
  gtf <- gtf_df %>% select(seqnames, strand, type, gene_id, transcript_id, exon_number, gene_name, ref_gene_id)
}

# Define functions to filter the df
novel <- function(gtf){
  gtf_novel <- gtf %>% filter(is.na(gene_name))
  return(gtf_novel)
}

annotated <- function(gtf){
  gtf_annotated <- gtf %>% filter(!is.na(gene_name))
  return(gtf_annotated)
}

transcripts <- function(gtf){
  gtf_transcripts <- gtf %>% filter(type == "transcript")
  return(gtf_transcripts)
}

exons <- function(gtf){
  gtf_exons <- gtf %>% filter(type == "exon")
  return(gtf_exons)
}

get_gene_ids <- function(gtf){
  gene_ids <- unique(gtf$gene_id)
  return(gene_ids)
}

get_gene_names <- function(gtf){
  gene_names <- unique(gtf$gene_name)
  return(gene_names)
}

get_transcript_ids <- function(gtf){
  transcript_ids <- unique(gtf$transcript_id)
    return(transcript_ids)
}

num_single_exons <- function(gtf){
  exon_one <- gtf %>% filter(exon_number == 1)
  exon_two <- gtf %>% filter(exon_number == 2)
  num_single_exon <- nrow(exon_one) - nrow(exon_two)
  return(num_single_exon)
}

num <- function(input){
  if(class(input) == "data.frame"){
    return(nrow(input))
  } else{
    return(length(input))
  }
}

############ Analyze Meta Assembly

gtf <- load_filter_gtf("../3_Assembly_2_MergeAssemblies/merged_assembly.gtf")

# Get number of transcripts in the meta-assembly
num(transcripts(gtf))
# 221417
num(transcripts(annotated(gtf)))
# 208311
num(transcripts(novel(gtf)))
# 13106

# Get number of exons in the meta-assembly
num(exons(gtf))
# 1367146
num(exons(annotated(gtf)))
# 1231077
num(exons(novel(gtf)))
# 136069

# Get number of transcripts composed of just a single exon
num_single_exons(gtf)
# 29390
num_single_exons(annotated(gtf))
# 28583
num_single_exons(novel(gtf))
# 807

# Get number of genes; gene_name comes from annotation, gene_id from assembly
num(get_gene_ids(gtf))
# 62375
num(get_gene_names(annotated(gtf)))
# 58646
num(get_gene_names(novel(gtf)))
# 1 = NA --> sanity check
num(get_gene_ids(annotated(gtf)))
# 61105
num(get_gene_ids(novel(gtf)))
# 6988
