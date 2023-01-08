#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --job-name="mapping"

module add UHTS/Aligner/hisat/2.2.1
# found through module avail UHTS

parent_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob

cell_line=$1

fastq_dir=${parent_dir}/fastq_files
sample_fwd=${fastq_dir}/${cell_line}*R1*.fastq.gz
sample_rev=${fastq_dir}/${cell_line}*R2*.fastq.gz

out_dir=${parent_dir}/2_Mapping/hisat2_out
out_file=${out_dir}/${cell_line}_hisat2.sam

ref_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/ref_genome/
ref_genome=${ref_dir}/*.fa
genome_name=GRCh38_genome

# build reference hisat2 index; $ref_genome points to the fasta file of the reference genome; $genome_name is the name we give the genome (used for generating the index)
# hisat2-build -p 16 $ref_genome $genome_name

# run hisat2 using the reference index just created; $sample_fwd and $sample_rev point to the fasta files to be mapped; $out_file is just the name we give the sam-file
ref_idx=${ref_dir}/${genome_name}
hisat2 -p 8 --dta -x $ref_idx -1 $sample_fwd -2 $sample_rev -S $out_file

# Options entered here are:
# ‘-p 8’ denoting the use of 8 threads
# ‘-–dta’ is used to generate output SAM files tailored for assemblers including StringTie
# ‘-x’ is used to denote the indexed reference genome
# ‘-1’ and ‘-2’ are used to denote our fwd and rev samples in a paired-end alignment
# ‘-S’ is used to denote that we would like our output in SAM format
# -u/--upto <int> align the first <int> reads or read pairs from the input (after the -s/--skip reads or pairs have been skipped), then stop. Default: no limit.