#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --job-name="NumReads"

fastq_dir=parent_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/fastq_files

#get number reads for the linked fastq files
for file in `ls -1 $fastq_dir/*.fastq.gz`;
do echo $file >> NumReads.txt;
zcat $file | grep -E '^@' | wc -l >> NumReads.txt;
done