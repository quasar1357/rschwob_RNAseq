#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=ProtCodPot
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/ProtCodPot_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/ProtCodPot_%j.e

# Add the modules
module load SequenceAnalysis/GenePrediction/cpat/1.2.4

course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
input=${course_dir}/results/6_IntAn_3_ProtCodPot/merged_assembly.bed
reference=${course_dir}/references/GRCh38.genome.fa
output=${course_dir}/results/6_IntAn_3_ProtCodPot/cpat_out.tab
hexamers=${course_dir}/results/6_IntAn_3_ProtCodPot/Human_Hexamer.tsv
logitModel=${course_dir}/results/6_IntAn_3_ProtCodPot/Human_logitModel.RData

cpat.py -g ${input} -r ${reference} -o ${output} -x ${hexamers} -d ${logitModel}