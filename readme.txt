Roman Schwob (roman.schwob@students.unibe.ch)

========================================================
=== Bioinformatics: Analysis of RNA sequencing reads ===
========================================================

This project is part of the course "RNA sequencing" (467713) of the University of Bern, taking place in Fall Semester 2022.

As part of group 1, I am analyzing the reads of the following cell lines:

--- --- --- --- --- --- --- ---  --- --- 

Holoclonal = 1.1, 1.2, 1.5

    Files:
   
    1_1_L3_R1_001_ij43KLkHk1vK.fastq.gz
    1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz
    
    1_2_L3_R1_001_DnNWKUYhfc9S.fastq.gz
    1_2_L3_R2_001_SNLaVsTQ6pwl.fastq.gz
    
    1_5_L3_R1_001_iXvvRzwmFxF3.fastq.gz
    1_5_L3_R2_001_iXCMrktKyEh0.fastq.gz

Parental = P1, P2, P3

    Files:
    
    P1_L3_R1_001_9L0tZ86sF4p8.fastq.gz
    P1_L3_R2_001_yd9NfV9WdvvL.fastq.gz
    
    P2_L3_R2_001_06FRMIIGwpH6.fastq.gz
    P2_L3_R2_001_06FRMIIGwpH6.fastq.gz
    
    P3_L3_R1_001_fjv6hlbFgCST.fastq.gz
    P3_L3_R2_001_xo7RBLLYYqeu.fastq.gz


--- --- --- --- --- --- --- ---  --- --- 


Overview:

1)  Read quality and statistics
    Goal:       How is the quality of the reads?
                How many reads do we have for each replicate?
    Software:   FASTQC V 0.11.9
                bash (grep)
    Scripts:    1_QC_1_FastQC.slurm
                1_QC_2_NumReads.slurm
    Input:      fastq files (see above)
    Output:     results of FastQC
                text file with number reads for confirmation

2)  Read mapping
    Goal:       Mapping reads onto human reference genome
                Alignment rates for the samples?
    Software:   HISAT2 V 2.2.1 (alternative = STAR)
                samtools V 1.10
    Reference:  Human genome version hg38/GRCh38
                Gencode, release 21, comprehensive gene annotation: https://www.gencodegenes.org/human/release_21.html
    Scripts:    2_Map_1_Index_RefGen.slurm
                2_Map_2_MapReads.slurm
                2_Map_3_SAM_index.slurm
                2_Map_4_BAM_to_SAM.slurm
    Input:      fastq files, forward and reverse each replicate
                reference genome (see above)
    Output:     BAM file for every replicate (sorted and indexed)

3)  Transcriptome assembly
    Goal:       How many exons, transcripts and genes are in your meta-assembly?
                How many of these are novel, i.e. do not have an associated GENCODE identifier?
                How many transcripts and genes are composed of just a single exon?
    Software:   StringTie (alternative = Scallop)
    Script:
    Input:      6 BAM files (1 of each cell line)
    Output:     One meta-assembly GTF format file (merged through stringtie --merge from 6 separate GTF files)

4)  Quantification
    Goal:       What units of expression are you using?
                Does the entire expression level across all genes add up to the expected amount?
                How many transcripts and genes did you detect?
                How many novel transcripts and genes did you detect?
    Software:   Kallisto (alternative = htseq-count)
    Script:
    Input: 
    Output:     Transcript and gene level expression tables

5)  Differential expression
    Goal:       Do known/expected genes change as expected?
    Software:   Sleuth (alternative = DESeq2)
    Script:
    Input: 
    Output:     Transcript- and gene-level differential expression tables

6)  Integrative analysis
    Goal:       How good are the 5’ and 3’ annotations of your transcripts?
                What percent of your novel transcripts are protein coding?
                How many novel “intergenic” genes have you identified?
    Software:   CPAT or CPCs
    Script:
    Input: 
    Output:     Statistics and plots addressing key questions

7)  Prioritization (Optional)
    Goal:       How would you prioritize your data to provide her with a ranked list of candidates?
    Software:
    Script:
    Input: 
    Output:     Ranked list of gene candidates