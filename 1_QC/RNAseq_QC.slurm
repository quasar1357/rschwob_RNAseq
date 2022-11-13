#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --job-name="QC"

module add UHTS/Quality_control/fastqc/0.11.9

parent_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob

#define and create result directory
res_dir=$parent_dir/QCres
mkdir $res_dir

#link fastq raw data in a separate folder; only necessary if not already linked
# fastq_dir=$parent_dir/fastq_files
# mkdir $fastq_dir
# ln -s /data/courses/rnaseq_course/lncRNAs/fastq/1* $fastq_dir
# ln -s /data/courses/rnaseq_course/lncRNAs/fastq/P* $fastq_dir

#run fastq on the linked fastq files
fastqc -t 2 -o $res_dir $fastq_dir/*.fastq.gz