#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=SAM2BAM
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/SAM2BAM_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/SAM2BAM_%j.e

# Add the modules
module add UHTS/Analysis/samtools/1.10

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_idx_dir=${course_dir}/results/2_Map_3_Index_RefGen_samtools
sam_dir=${course_dir}/results/2_Map_2_MapReads
output_dir=${course_dir}/results/2_Map_4_SAM2BAM
    mkdir ${output_dir}

# Define the files to be used and created
samtools_idx=${ref_idx_dir}/*.fai
    # The samtools index created in the previous step
cell_line=$1
    # The cell line is given as a parameter to the script
sam_file=${sam_dir}/${cell_line}*.sam
    # Use the SAM file belonging to the cell line
bam_unsorted=${output_dir}/${cell_line}_unsorted.bam
bam_sorted=${output_dir}/${cell_line}_sorted.bam
    # The output BAM files

# Convert SAM to BAM
samtools view -b -t ${samtools_idx} ${sam_file} > ${bam_unsorted}
    # Options:
        # "-b": output as BAM
        # "-t": samtool index (FILE listing reference names and lengths)

# Sort BAM
samtools sort -o ${bam_sorted} ${bam_unsorted}
    # Options:
        # "sort": Sort alignments by leftmost coordinates, or by read name when -n is used. 
        # "-o":Write final output to FILE rather than standard output
# rm ${bam_unsorted}

# Index BAM
samtools index ${bam_sorted}
