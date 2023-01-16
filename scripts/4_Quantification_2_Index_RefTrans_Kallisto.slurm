#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=Index_RefTrans_Kallisto
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Index_RefTrans_Kallisto_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Index_RefTrans_Kallisto_%j.e

# Add the modules
module add UHTS/Analysis/kallisto/0.46.0

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
reftrans_in_dir=${course_dir}/results/4_Quantification_1_Create_RefTrans
output_dir=${course_dir}/results/4_Quantification_2_Index_RefTrans_Kallisto
    mkdir ${output_dir}

# Define the files to be used and created
reftrans_in=${reftrans_in_dir}/ref_transcriptome.fasta
    # Define the reference transcriptome file for input
kallisto_index=${output_dir}/ref_trans_kallisto_index.idx
    # What file shall the Kallisto index be written to?

# Create kallisto index
kallisto index -i ${kallisto_index} ${reftrans_in}
    # Options:
        # "index": Makes an index for Kallisto analysis
        # "-i": Filename for the kallisto index to be constructed
        # "${reftrans_in}": fasta file to be indexed (i.e. reference transcriptome)