Roman Schwob (roman.schwob@students.unibe.ch)

========================================================
=== Bioinformatics: Analysis of RNA sequencing reads ===
========================================================

This project is part of the course "RNA sequencing" (467713) of the University of Bern, taking place in Fall Semester 2022/2023.

As part of subgroup 1 in the lncRNA group, I analyzed the reads of the holoclonal and compared them to the parental cell lines.


--- --- --- --- --- --- --- --- --- ---
Datasets
--- --- --- --- --- --- --- --- --- ---

Holoclonal (3 replicates: 1_1, 1_2, 1_5), reads from TruSeq stranded libraries

    Files (3 replicates, each with R1 and R2 = reverse and forward, respectively):
   
    1_1_L3_R1_001_ij43KLkHk1vK.fastq.gz
    1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz
    
    1_2_L3_R1_001_DnNWKUYhfc9S.fastq.gz
    1_2_L3_R2_001_SNLaVsTQ6pwl.fastq.gz
    
    1_5_L3_R1_001_iXvvRzwmFxF3.fastq.gz
    1_5_L3_R2_001_iXCMrktKyEh0.fastq.gz

Parental (3 replicates: P1, P2, P3), reads from TruSeq stranded libraries

    Files (3 replicates, each with R1 and R2 = reverse and forward, respectively):
    
    P1_L3_R1_001_9L0tZ86sF4p8.fastq.gz
    P1_L3_R2_001_yd9NfV9WdvvL.fastq.gz
    
    P2_L3_R2_001_06FRMIIGwpH6.fastq.gz
    P2_L3_R2_001_06FRMIIGwpH6.fastq.gz
    
    P3_L3_R1_001_fjv6hlbFgCST.fastq.gz
    P3_L3_R2_001_xo7RBLLYYqeu.fastq.gz

Reference genome

    Human genome, version GRCh38
    GENCODE release 21 (https://www.gencodegenes.org/human/release_21.html):
    - Comprehensive gene annotation with ALL regions in GTF format
    - Genome sequence (GRCh38) with ALL regions in fasta format
    Uni Basel, Human Build 38 (https://www.polyasite.unibas.ch/atlas):
    - Atlas BED file with PolyA sites 
    FANTOM5 collection (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/):
    - hg38_fair+new_CAGE_peaks_phase1and2 BED file
    

--- --- --- --- --- --- --- --- --- ---
Data analysis steps
--- --- --- --- --- --- --- --- --- ---

1)  Read quality and statistics
    Goal:       How is the quality of the reads?
                How many reads do we have for each replicate?
    Software:   FastQC 0.11.9
                bash (grep)
    Scripts:    1_QC_1_FastQC.slurm
                1_QC_2_NumReads.slurm
    Input:      fastq files
    Output:     Results of FastQC
                Text file with number reads for confirmation

2)  Read mapping
    Goal:       Mapping reads onto human reference genome
                Alignment rates for the samples?
    Software:   HISAT2 2.2.1 (alternative = STAR)
                Samtools 1.10
                IGV 2.15.4
    Scripts:    2_Map_1_Index_RefGen_hisat2.slurm
                2_Map_2_MapReads.slurm
                2_Map_3_Index_RefGen_samtools.slurm
                2_Map_4_SAM2BAM.slurm
    Input:      fastq files, reverse and forward each replicate
                Reference genome in fasta format
    Output:     BAM file for every replicate (sorted and indexed)
                Text file with the alignment rates (in 2_Map_2_MapReads)

3)  Transcriptome assembly
    Goal:       How many exons, transcripts and genes are in the meta-assembly?
                How many transcripts are composed of just a single exon?
                How many of these are novel, i.e. do not have an associated GENCODE identifier?
    Software:   StringTie 1.3.3b (alternative = Scallop)
                R 4.2.2
    Scripts:    3_Assembly_1_SingleAssembly.slurm
                3_Assembly_2_MergeAssemblies.slurm
                3_Assembly_3_FilterCount_MetaAssembly.R
    Input:      6 BAM files (1 of each cell line)
                Reference genome in gtf format
    Output:     One meta-assembly GTF format file (merged through stringtie --merge from 6 separate GTF files)
                Txt file with numbers of genes, transcripts, exons that are novel, annotated, single-exon

4)  Quantification
    Goal:       What units of expression are we using?
                Does the entire expression level across all genes add up to the expected amount?
    Software:   Cufflinks 2.2.1 (to generate reference transcriptome)
                kallisto 0.46.0 (alternative = htseq-count)
                R 4.2.2
    Scripts:    4_Quantification_1_Create_RefTrans.slurm
                4_Quantification_2_Index_RefTrans_Kallisto.slurm
                4_Quantification_3_ExprLevels.slurm
                4_Quantification_4_Validation.R
    Input:      Meta-assembly GTF
                Reference genome in fasta format
                fastq files, reverse and forward each replicate
    Output:     Absolute expression tables (abundance.h5)
                Txt file with total tpm (= 1'000'000) for validation of expression levels
                Txt file with the number of expressed transcripts for each replicate

5)  Differential expression
    Goal:       Do known/expected genes change in expression as expected?
    Software:   sleuth 0.30.1 (alternative = DESeq2)
                R 4.2.2
    Scripts:    5_DiffExpr_2_Gene_Transcript_Map.R
                5_DiffExpr_3_DiffExpr_TransLevel.R
                5_DiffExpr_4_DiffExpr_GeneLevel.R
    Input:      Absolute expression tables (abundance.h5 from kallisto)
                (Meta-assembly GTF for gene-transcript map)
    Output:     Transcript and gene level differential expression tables & plots

6)  Integrative analysis
    Goal:       How good are the 5’ and 3’ annotations of your transcripts?
                How many novel “intergenic” genes have you identified?
                What percent of your novel transcripts are protein coding?
    Software:   R 4.2.2
                BEDTools 2.29.2
                CPAT 1.2.4 and CPC 2.0
    Scripts:    6_IntAn_1_CreateBed.R
                6_IntAn_2_Find_TSS_PolyA_intergenic.slurm
                6_IntAn_3_ProtCodPot.slurm
    Input:      Merged meta-assembly GTF
                BED references for polyA and TSS
                Reference genome in gtf format (for finding "intergenic" genes)
                Human hexamer frequencies and logitModel for CPAT (provided by sourceforge)
    Output:     Statistics and plots addressing key questions

7)  Summary/prioritization (Optional)
    Goal:       How would you prioritize your data to provide it with a ranked list of candidates?
    Software:   R 4.2.2
    Script:     7_1_Create_Summary.R
                7_2_Analyse_Results.R
    Input:      Merged meta-assembly GTF
                DiffExpr tables
                TSS, polyA and intergenic tables
                Protein coding potential table
                (Reference genome in gtf format)
    Output:     Ranked list of gene candidates


--- --- --- --- --- --- --- --- --- ---
Software used
--- --- --- --- --- --- --- --- --- ---

FastQC:		    https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
HISAT2:		    https://daehwankimlab.github.io/hisat2/manual/
Samtools:	    https://www.htslib.org/
IGV:		    https://igv.org/; https://software.broadinstitute.org/software/igv/
StringTie:	    https://ccb.jhu.edu/software/stringtie/
Cufflinks:	    https://cole-trapnell-lab.github.io/cufflinks/
kallisto:	    https://pachterlab.github.io/kallisto/about.html
BEDTools:       https://bedtools.readthedocs.io/en/latest/
CPAT 1.2.4:     https://rna-cpat.sourceforge.net/
(CPAT 3         https://cpat.readthedocs.io/en/latest/ )
(CPC 2.0:        https://cpc2.gao-lab.org/ )

rtracklayer:	https://rdrr.io/bioc/rtracklayer/
sleuth:		    https://pachterlab.github.io/sleuth/about

Overview over all kinds of bioinformatics software, including most of the above:
https://bioinformaticshome.com/
Overview over software on the ibu cluster:
https://www.vital-it.ch/
