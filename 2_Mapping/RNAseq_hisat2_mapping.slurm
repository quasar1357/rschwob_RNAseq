#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --job-name="mapping"

module add UHTS/Aligner/hisat/2.2.1
module add UHTS/Analysis/samtools/1.10

# Define and create paths and directories
parent_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob

cell_line=$1

fastq_dir=${parent_dir}/fastq_files
sample_fwd=${fastq_dir}/${cell_line}*R1*.fastq.gz
sample_rev=${fastq_dir}/${cell_line}*R2*.fastq.gz

hisat2_out_dir=${parent_dir}/2_Mapping/hisat2_out
sam_file=${hisat2_out_dir}/${cell_line}_hisat2.sam

ref_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/ref_genome/
ref_genome=${ref_dir}/*.fa
genome_name=GRCh38_genome

# build reference hisat2 index; $ref_genome points to the fasta file of the reference genome; $genome_name is the name we give the genome (used for generating the index)
####### ####### ####### ####### ####### ####### 
# hisat2-build -p 16 $ref_genome $genome_name
####### ####### ####### ####### ####### #######  

# run hisat2 using the reference index just created; $sample_fwd and $sample_rev point to the fasta files to be mapped; $sam_file is just the name we give the sam-file
####### ####### ####### ####### ####### ####### 
# hisat2_idx=${ref_dir}/${genome_name}
# hisat2 -p 8 --dta -x $hisat2_idx -1 $sample_fwd -2 $sample_rev -S $sam_file
####### ####### ####### ####### ####### ####### 

# Options entered here are:
# ‘-p 8’ denoting the use of 8 threads
# ‘-–dta’ is used to generate output SAM files tailored for assemblers including StringTie
# ‘-x’ is used to denote the indexed reference genome
# ‘-1’ and ‘-2’ are used to denote our fwd and rev samples in a paired-end alignment
# ‘-S’ is used to denote that we would like our output in SAM format
# -u/--upto <int> align the first <int> reads or read pairs from the input (after the -s/--skip reads or pairs have been skipped), then stop. Default: no limit.


### CONVERSION TO BAM ###

# Define and create paths and direcories
samtools_out_dir=${parent_dir}/2_Mapping/samtools_out
mkdir ${samtools_out_dir}
sam_idx=${samtools_out}/${genome_name}_samtools_idx.fai

bam_unsorted=${samtools_out_dir}/${cell_line}_unsorted.bam
bam_sorted=${samtools_out_dir}/${cell_line}_sorted.bam

# Index reference for samtools (only necessasry to do once)
####### ####### ####### ####### ####### ####### 
# samtools faidx ${ref_genome} > ${sam_idx}
####### ####### ####### ####### ####### ####### 

# Convert SAM to BAM
####### ####### ####### ####### ####### ####### 
samtools view -b -t ${sam_idx} ${sam_file} > ${bam_unsorted}
####### ####### ####### ####### ####### ####### 
# Options entered here are:
    #"-b": output BAM
    #"-t": FILE listing reference names and lengths

# Sort BAM
####### ####### ####### ####### ####### ####### 
samtools sort -o ${bam_sorted} ${bam_unsorted}
####### ####### ####### ####### ####### ####### 
# Options entered here are:
    #"-o":Write final output to FILE rather than standard output

# Index BAM
####### ####### ####### ####### ####### ####### 
samtools index ${bam_sorted}
####### ####### ####### ####### ####### ####### 