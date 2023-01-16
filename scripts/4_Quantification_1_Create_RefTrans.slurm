#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=Find_TSS_PolyA_intergenic
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Find_TSS_PolyA_intergenic_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Find_TSS_PolyA_intergenic_%j.e

# Add the modules
module add UHTS/Assembler/cufflinks/2.2.1

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_dir=${course_dir}/references
gtf_in_dir=${course_dir}/results/3_Assembly_2_MergeAssemblies
output_dir=${course_dir}/results/4_Quantification_1_Create_RefTrans
    mkdir ${output_dir}
    cd ${output_dir}

# Define the files to be used and created
gtf_in=${gtf_in_dir}/merged_assembly.gtf
    # Define the gtf file for input
ref_fasta=${ref_dir}/GRCh38.genome.fa
    # Where is the reference fasta file of the reference genome?
    # NOTE: cufflinks will build a fasta index automatically if there is none in the directory

# Make the transcriptome fasta with gffread for the kallisto index
gffread ${gtf_in} -g ${ref_fasta} -w ref_transcriptome.fasta
    # Options:
        # "${gtf_in}": Input GTF file (merged meta-assembly).
        # "-g": full path to a multi-fasta file with the genomic sequences for all input mappings
        # "-w": write a fasta file with spliced exons for each GFF transcript