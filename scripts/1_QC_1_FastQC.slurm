#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=FastQC
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/FastQC_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/FastQC_%j.e

#Add the modules
module add UHTS/Quality_control/fastqc/0.11.9

#Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
fastq_dir=${course_dir}/fastq_files
output_dir=${course_dir}/script_outputs/1_QC_1_FastQC
    mkdir ${output_dir}

#Make the quality analysis
fastqc -t 2 -o ${output_dir} ${fastq_dir}/*.fastq.gz
    #Options entered here are:
        #"-t 2": Specifies the number of files which can be processed simultaneously.
        #"-o ${output_dir}": Create all output files in the specified output directory.
        #"${fastq_dir}/*.fastq.gz": Input all fastq files in the input directory.