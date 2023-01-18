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

# Define directories and files to use or create
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
# input=${course_dir}/results/6_IntAn_3_ProtCodPot/merged_assembly_filtered.bed
# input=${course_dir}/results/6_IntAn_3_ProtCodPot/merged_assembly.bed
input=${course_dir}/results/4_Quantification_1_Create_RefTrans/ref_transcriptome.fasta
reference=${course_dir}/references/GRCh38.genome.fa
output=${course_dir}/results/6_IntAn_3_ProtCodPot/cpat_out
hexamers=${course_dir}/results/6_IntAn_3_ProtCodPot/Human_Hexamer.tsv
logitModel=${course_dir}/results/6_IntAn_3_ProtCodPot/Human_logitModel.RData

# Run CPAT
cpat.py -g ${input} -o ${output} -x ${hexamers} -d ${logitModel}
# Options:
    # -g GENE_FILE: Transcripts either in BED format or mRNA sequences in FASTA format: If this is BED format file, ‘-r’ must be specified; if this is mRNA sequence file in FASTA format, ignore the ‘-r’ option. The input BED or FASTA file could be regular text file or compressed file (*.gz, *.bz2) or accessible url.
    # -o OUT_FILE: output file. Tab separated text file: geneID <tab> mRNA size <tab> ORF size <tab> Fickett Score <tab> Hexamer Score<tab>Coding Probability.
    # -x HEXAMER_DAT: Prebuilt hexamer frequency table (Human, Mouse, Fly, Zebrafish). Run ‘make_hexamer_tab.py’ to make this table out of your own training dataset.
    # -d LOGIT_MODEL: Prebuilt training model (Human, Mouse, Fly, Zebrafish). Run ‘make_logitModel.py’ to build logit model out of your own training datset

