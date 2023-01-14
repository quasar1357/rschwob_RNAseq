library(tidyverse)
library(rtracklayer)
library(sleuth)



# Load prepared, required files
experiment_info <- read.csv("../5_DiffExpr_1_Experiment_Table/experiment_table.csv")
gene_transcript_map <- read.csv(file = "../5_DiffExpr_2_Gene_Transcript_Map/gene_transcript_map.csv")


########## Construct "sleuth object", store information about experiment and details of model

### 1.Load kallisto processed data into the object
so_gene <- sleuth_prep(experiment_info, transformation_function = function(x) log2(x + 0.5), target_mapping = gene_transcript_map, aggregation_column = "gene_id", gene_mode = TRUE, read_bootstrap_tpm = TRUE , extra_bootstrap_summary = TRUE)
# Output:
  # 78930 targets passed the filter
  # 23518 genes passed the filter
  # 20 warnings (see txt file for details)
# Options:
  # "transformation_function": By default the transformation of counts is natural log, which would make the output fold changes somewhat more difficult to interpret. By specifying the transformation_function to be log2(x + 0.5) we are ensuring our output fold changes are log2.
  # "target_mapping": a data.frame that has at least one column 'target_id' and others that denote the mapping for each target.
    # If it is not NULL, target_mapping is joined with many outputs where it might be useful. For example, you might have columns 'target_id', 'ensembl_gene' and 'entrez_gene' to denote different transcript to gene mappings.
    # Note that sleuth_prep will treat all columns as having the 'character' data type.
  # "aggregation_column": a string of the column name in target_mapping to aggregate targets (typically to summarize the data on the gene level).
    # The aggregation is done using a p-value aggregation method when generating the results table. See sleuth_results for more information.
  # "gene_mode": Set this to TRUE to get the old counts-aggregation method for doing gene-level analysis. This requires aggregation_column to be set.
    # If TRUE, this will override the p-value aggregation mode, but will allow for gene-centric modeling, plotting, and results.
  # "read_bootstrap_tpm": read and compute summary statistics on bootstraps on the TPM. This is not necessary for typical analyses; it is only needed for some plots (e.g. plot_bootstrap) and if TPM values are used for sleuth_fit. Default is FALSE.
  # "extra_bootstrap_summary": if TRUE, compute extra summary statistics for estimated counts. This is not necessary for typical analyses; it is only needed for certain plots (e.g. plot_bootstrap). Default is FALSE.


### 2. Estimate parameters for sleuth response error measurement (full) model
so_gene <- sleuth_fit(so_gene, ~condition, "full")
# Output:
  # 5 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
  # The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
  # These are the target ids with NA values: AC016907.3, LRP2, RNF152, RP13-580B18.4, ZDHHC22  

### 3. Estimate parameters for sleuth reduced model
so_gene <- sleuth_fit(so_gene, ~1, "reduced")
# Output:
  # 6 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
  # The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
  # These are the target ids with NA values: AC016907.3, LRP2, RNF152, RP13-580B18.4, ZDHHC22, MT-CO1

### The models that have been fit can always be examined with the models() function:
# models(so_transcript)


########## Perform different types of analyses on the sleuth object

### 1. Perform differential analysis (testing) using the likelihood ratio test
so_gene <- sleuth_lrt(so_gene, "reduced", "full")

### 2. Perform differential analysis (testing) using the Wald test
so_gene <- sleuth_wt(so_gene, "conditionparental")


########## Examine results of test

sleuth_lrt_gene <- sleuth_results(so_gene, "reduced:full", "lrt", show_all = FALSE)
sleuth_lrt_gene_significant <- filter(sleuth_lrt_gene, qval <= 0.05)
sleuth_wt_gene <- sleuth_results(so_gene, "conditionparental", show_all = TRUE)
sleuth_wt_gene_significant <- filter(sleuth_wt_gene, qval <= 0.05)


########## Save into files

out_dir <- "../5_DiffExpr_4_DifExpr_GeneLevel"
write.csv(sleuth_lrt_gene, file = paste0(out_dir, "/sleuth_lrt_gene.csv"), row.names = FALSE)
write.csv(sleuth_lrt_gene_significant, file = paste0(out_dir, "/sleuth_lrt_gene_significant.csv"), row.names = FALSE)
write.csv(sleuth_wt_gene, file = paste0(out_dir, "/sleuth_wt_gene.csv"), row.names = FALSE)
write.csv(sleuth_wt_gene_significant, file = paste0(out_dir, "/sleuth_wt_gene_significant.csv"), row.names = FALSE)
sleuth_save(so_gene, file = paste0(out_dir, "/so_gene"))


########## Create plots

# Read in the sleuth object (only necessary if not done in same session and don't want to run everything again)
in_dir <- "../5_DiffExpr_4_DifExpr_GeneLevel"
so_gene <- sleuth_load(paste0(in_dir, "/so_gene"))

# Show sleuth live (tool for plots)
sleuth_live(so_gene)


########## Read in the raw data (only needed to analyse if don't want to rerun the whole script)

# in_dir <- "../5_DiffExpr_4_DifExpr_GeneLevel"
# sleuth_lrt_gene <- read.csv(file = paste0(in_dir, "/sleuth_lrt_gene.csv"))
# sleuth_lrt_gene_significant <- read.csv(file = paste0(in_dir, "/sleuth_lrt_gene_significant.csv"))
# sleuth_wt_gene <- read.csv(file = paste0(in_dir, "/sleuth_wt_gene.csv"))
# sleuth_wt_gene_significant <- read.csv(file = paste0(in_dir, "/sleuth_wt_gene_significant.csv"))