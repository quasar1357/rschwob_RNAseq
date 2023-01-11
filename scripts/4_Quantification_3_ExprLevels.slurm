#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=ExprLevels
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/ExprLevels_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/ExprLevels_%j.e

# Add the modules
module add UHTS/Analysis/kallisto/0.46.0

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
fastq_dir=${course_dir}/fastq_files
reftrans_idx_dir=${course_dir}/results/4_Quantification_2_Index_RefTrans_Kallisto
merged_assembly_dir=${course_dir}/results/3_Assembly_2_MergeAssemblies
output_dir=${course_dir}/results/4_Quantification_3_ExprLevels
    mkdir ${output_dir}
cell_line=$1
    # The cell line is given as a parameter to the script
cell_line_out_dir=${output_dir}/${cell_line}
    mkdir ${cell_line_out_dir}

# Define the files to be used and created
sample_fwd=${fastq_dir}/${cell_line}*R1*.fastq.gz
sample_rev=${fastq_dir}/${cell_line}*R2*.fastq.gz
    # fastq files with the forward and reverse reads (original read fastq files of the sample cell lines)
reftrans_idx=${reftrans_idx_dir}/ref_trans_kallisto_index.idx
    # Define the reference transcriptome index file for input
merged_assembly=${merged_assembly_dir}/merged_assembly.gtf
    # What file shall the Kallisto index be written to?


#Make the kallisto quantification
    kallisto quant -i ${reftrans_idx} -o ${cell_line_out_dir} --rf-stranded -b 30 -g ${merged_assembly} ${sample_fwd} ${sample_rev}
    # Options:
        # "quant": runs the quantification algorithm.
        # "–i": Filename for the indexed reference transcriptome (kallisto index) to be used for quantification
        # "–o": Directory to write output to.
        # "--rf-stranded": Strand specific reads, first read reverse.
        # "-b 30": Number of bootstrap samples
        # "–g": GTF file for transcriptome information (merged assembly gtf)
        # "${sample_fwd}/${sample_rev}": file (fastq) of forward and backward reads.