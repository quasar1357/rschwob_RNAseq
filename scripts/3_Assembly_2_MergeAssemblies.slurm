#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --time=05:00:00
#SBATCH --job-name=MergeAssemblies
#SBATCH --mail-user=roman.schwob@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --output=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/MergeAssemblies_%j.o
#SBATCH --error=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob/slurm_out_files/MergeAssemblies_%j.e

# Add the modules
module add UHTS/Aligner/stringtie/1.3.3b

# Define (and create) input and output directories
course_dir=/data/courses/rnaseq_course/lncRNAs/Project1/users/rschwob
ref_dir=${course_dir}/references
gtf_in_dir=${course_dir}/results/3_Assembly_1_SingleAssembly
output_dir=${course_dir}/results/3_Assembly_2_MergeAssemblies
    mkdir ${output_dir}

# Define the files to be used and created
all_gtf_in=${gtf_in_dir}/*.gtf
    # Define all gtf single assembly files for input
gtf_list=${output_dir}/single_GTF_list.txt
    # A list of input GTF files will be created; what file should this be written to?
ref_gtf=${ref_dir}/gencode.v21.chr_patch_hapl_scaff.annotation.gtf
    # Where is the reference GTF file
gtf_out=${output_dir}/merged_assembly.gtf
    # What gtf file should be written to

# Merge single assemblies into one meta-assembly GTF
ls -1 ${all_gtf_in} > ${gtf_list}
    # Create a list of all input GTF files
stringtie --merge --rf -o ${gtf_out} -G ${ref_gtf} ${gtf_list}
    # Do the actual merge of the assemblies
    # Options:
        # "--merge": Transcript merge mode. This is a special usage mode of StringTie, distinct from the assembly usage mode described above.
            # In the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts.
            # This mode is used in the new differential analysis pipeline to generate a global, unified set of transcripts (isoforms) across multiple RNA-Seq samples.
        # "–-rf": Assumes a stranded library fr-firststrand.
        # "-o: output file name for the merged transcripts GTF.
        # "-G: reference annotation to include in the merging.
        # ${gtf_list}: text file with a list (one per line) of GTF files that you would like to merge together into a single GTF file.