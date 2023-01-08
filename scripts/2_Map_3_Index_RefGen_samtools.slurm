#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=Index_RefGen_samtools
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Index_RefGen_samtools_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Index_RefGen_samtools_%j.e

# Add the modules
module add UHTS/Analysis/samtools/1.10

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_dir=${course_dir}/references
output_dir=${course_dir}/results/2_Map_3_Index_RefGen_samtools
    mkdir ${output_dir}

# Index reference genome with/for samtools
ref_fasta=${ref_dir}/*.fa
    # $ref_fasta points to the fasta file of the reference genome (should in our case just be one fasta file of GRCh38_genome)
genome_name=GRCh38_genome
sam_index=${output_dir}/${genome_name}_samtools_idx.fai
    # $sam_index is the output name for the samtools index of the reference genome
samtools faidx ${ref_fasta} > ${sam_index}
