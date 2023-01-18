#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=Find_TSS_PolyA_intergenic
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Find_TSS_PolyA_intergenic_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/Find_TSS_PolyA_intergenic_%j.e

# Add the modules
module add UHTS/Analysis/BEDTools/2.29.2
module add SequenceAnalysis/GenePrediction/cpat/1.2.4;

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_dir=${course_dir}/references
created_bed_dir=${course_dir}/results/6_IntAn_1_CreateBed
output_dir=${course_dir}/results/6_IntAn_2_Find_TSS_PolyA_intergenic
    mkdir ${output_dir}

# Define the files to be used and created
ref_TSS_bed=${ref_dir}/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
ref_polyA_bed=${created_bed_dir}/atlas.clusters.2.0.GRCh38.96_ADAPTED.bed
ref_total_bed=${created_bed_dir}/reference_TOTAL.bed

assembly_TSS_bed=${created_bed_dir}/merged_assembly_TSS.bed
assembly_polyA_bed=${created_bed_dir}/merged_assembly_polyA.bed
assembly_total_bed=${created_bed_dir}/merged_assembly_TOTAL.bed

output_TSS=${output_dir}/TSS_intersect.tsv
output_polyA=${output_dir}/polyA_intersect.tsv
output_total=${output_dir}/total_intersect.tsv

# Find intersect of the TSS
bedtools intersect -a ${assembly_TSS_bed} -b ${ref_TSS_bed} -u -s > ${output_TSS}
    # Options:
        # "-a": BAM/BED/GFF/VCF file "A". Each feature in A is compared to B in search of overlaps.
        # "-b": One or more BAM/BED/GFF/VCF file(s) "B".
        # "-s": Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
        # "-u": Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B.

# Find intersect of the polyA
bedtools intersect -a ${assembly_polyA_bed} -b ${ref_polyA_bed} -u -s > ${output_polyA}
    # Options:
        # "-a": BAM/BED/GFF/VCF file "A". Each feature in A is compared to B in search of overlaps.
        # "-b": One or more BAM/BED/GFF/VCF file(s) "B".
        # "-s": Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
        # "-u": Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B.

# Find "intergenic" entries by finding entries that do NOT OVERLAP with any entry in the reference gtf
bedtools intersect -a ${assembly_total_bed} -b ${ref_total_bed} -v -s > ${output_total}
    # Options:
        # "-a": BAM/BED/GFF/VCF file "A". Each feature in A is compared to B in search of overlaps.
        # "-b": One or more BAM/BED/GFF/VCF file(s) "B".
        # "-s": Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
        # "-v": Only report those entries in A that have no overlap in B.
