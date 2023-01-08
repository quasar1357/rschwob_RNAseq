#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=FastQC
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/NumReads_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/NumReads_%j.e

#Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
fastq_dir=${course_dir}/fastq_files
output_dir=${course_dir}/results/1_QC_2_NumReads
    mkdir ${output_dir}

#Get the number of reads for the different files
input_files=$(ls -1 $fastq_dir/*.fastq.gz)
for file in $input_files
    do
        nr_reads=$(zcat ${file} | grep -E '^@' | wc -l)
        echo "$(basename ${file}): ${nr_reads} reads" >> ${output_dir}/NumReads.txt
    done
