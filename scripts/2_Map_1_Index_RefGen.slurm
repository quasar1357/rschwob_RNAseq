#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=Index_RefGen
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Index_RefGen_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Index_RefGen_%j.e

# Add the modules
module add UHTS/Aligner/hisat/2.2.1

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_dir=${course_dir}/references
output_dir=${course_dir}/results/2_Map_1_Index_RefGen
    mkdir ${output_dir}
    cd ${output_dir}

# Build reference hisat2 index
ref_fasta=${ref_dir}/*.fa
    # $ref_fasta points to the fasta file of the reference genome (should in our case just be one fasta file of GRCh38_genome)
genome_name=GRCh38_genome
    # $genome_name is the name we give the genome (used for generating the index)
hisat2-build -p 16 $ref_fasta $genome_name
    # Options entered here are:
        #"-p 16": Launch 16 parallel build threads.
        #"${ref_fasta}": A comma-separated list of FASTA files containing the reference sequences to be aligned to.
        #"${genome_name}": The basename of the index files to write.
