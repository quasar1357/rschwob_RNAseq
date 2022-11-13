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


General steps:

1) Read quality and statistics
    Input: fastq files of reads (see above)
    Outputs: results of fastQC (folder QCres), text file with number reads for confirmation

2) Read mapping
    Input: 
    Output: BAM file for every replicate

3) Transcriptome assembly
    Input: 
    Output: One meta-assembly GTF format file

4) Quantification
    Input: 
    Output: Transcript and gene level expression tables

5) Differential expression
    Input: 
    Output: Transcript- and gene-level differential expression tables

6) Integrative analysis
    Input: 
    Output: Statistics and plots addressing key questions

7) Prioritization (Optional)
    Input: 
    Output: Ranked list of gene candidates