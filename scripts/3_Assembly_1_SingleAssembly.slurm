#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=SingleAssembly
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/SingleAssembly_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/SingleAssembly_%j.e

# Add the modules
module add UHTS/Aligner/stringtie/1.3.3b

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_dir=${course_dir}/references
bam_dir=${course_dir}/results/2_Map_4_SAM2BAM
output_dir=${course_dir}/results/3_Assembly_1_SingleAssembly
    mkdir ${output_dir}

# Define the files to be used and created
cell_line=$1
    # The cell line is given as a parameter to the script
bam_in=${bam_dir}/${cell_line}_sorted.bam
    # Use the BAM file belonging to the cell line
ref_gtf=${ref_dir}/gencode.v21.chr_patch_hapl_scaff.annotation.gtf
    # Where is the reference GTF file
gtf_out=${output_dir}/${cell_line}.gtf
    # What gtf file should be written to
tab_out=${output_dir}/${cell_line}_gene_abund.tab
    # What gene abundance tab file should be written to (see below, option -A)

# Make assembly for this cell line
stringtie --rf -G ${ref_gtf} -o ${gtf_out} ${bam_in} -A ${tab_out}
    # Options:
        # "--rf": Assumes a stranded library fr-firststrand.
        # "-G: Use a reference annotation file (in GTF or GFF3 format) to guide the assembly process. The output will include expressed reference transcripts as well as any novel transcripts that are assembled.
        # "-o: Sets the name of the output GTF file where StringTie will write the assembled transcripts.
        # ${bam_in}: Input BAM file.
        # "-A: Gene abundances will be reported (tab delimited format) in the output file with the given name.