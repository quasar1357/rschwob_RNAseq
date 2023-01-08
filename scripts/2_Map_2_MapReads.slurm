#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=MapReads
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/MapReads_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/MapReads_%j.e

# Add the modules
module add UHTS/Aligner/hisat/2.2.1

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_idx_dir=${course_dir}/results/2_Map_1_Index_RefGen_hisat2
fastq_dir=${course_dir}/fastq_files
output_dir=${course_dir}/results/2_Map_2_MapReads
    mkdir ${output_dir}

# Define the files to be used and created
genome_name=GRCh38_genome
    # $genome_name is the name of the genome (same as used for generating the index)
ref_idx=${ref_idx_dir}/${genome_name}
    # $ref_idx means the basename of the index files. By default, the index files are named NAME.1.ht2 where we need to give the basename "NAME" (= $ref_idx)
cell_line=$1
    # The cell line is given as a parameter to the script
sample_fwd=${fastq_dir}/${cell_line}*R1*.fastq.gz
sample_rev=${fastq_dir}/${cell_line}*R2*.fastq.gz
    # $sample_fwd and $sample_rev point to the fasta files to be mapped
sam_out_file=${output_dir}/${cell_line}_hisat2.sam
    # $sam_out_file is just the name (incl. path) we give the sam-file to be written

# Run hisat2 using the reference index created before
hisat2 -p 8 --dta --rna-strandness RF -x ${ref_idx} -1 ${sample_fwd} -2 ${sample_rev} -S ${sam_out_file}
    # Options entered here are:
        # ‘-p 8’ denoting the use of 8 threads
        # ‘-–dta’ is used to generate output SAM files tailored for assemblers including StringTie
        # ‘--rna-strandness RF’ specifies strand-specific information.
        # ‘-x’ is used to denote the indexed reference genome
        # ‘-1’ and ‘-2’ are used to denote our fwd and rev samples in a paired-end alignment
        # ‘-S’ is used to denote that we would like our output in SAM format