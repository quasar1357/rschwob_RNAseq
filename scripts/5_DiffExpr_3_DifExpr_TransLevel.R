library(tidyverse)
library(rtracklayer)
library(sleuth)


# Load prepared, required files
experiment_info <- read.csv("../results/5_DiffExpr_1_Experiment_Table/experiment_table.csv")


########## Construct "sleuth object", store information about experiment and details of model

### 1.Load kallisto processed data into the object
so_transcript <- sleuth_prep(experiment_info, transformation_function = function(x) log2(x + 0.5), read_bootstrap_tpm = TRUE , extra_bootstrap_summary = TRUE)
# Output:
  # 78930 targets passed the filter
# Options:
  # "transformation_function": By default the transformation of counts is natural log, which would make the output fold changes somewhat more difficult to interpret. By specifying the transformation_function to be log2(x + 0.5) we are ensuring our output fold changes are log2.
  # "read_bootstrap_tpm": read and compute summary statistics on bootstraps on the TPM. This is not necessary for typical analyses; it is only needed for some plots (e.g. plot_bootstrap) and if TPM values are used for sleuth_fit. Default is FALSE.
  # "extra_bootstrap_summary": if TRUE, compute extra summary statistics for estimated counts. This is not necessary for typical analyses; it is only needed for certain plots (e.g. plot_bootstrap). Default is FALSE.


### 2. Estimate parameters for sleuth response error measurement (full) model
so_transcript <- sleuth_fit(so_transcript, ~condition, "full")
# Output:
  # 7 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
  # The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
  # These are the target ids with NA values: ENST00000378610.1, ENST00000411386.1, ENST00000498450.1, ENST00000534689.3, ENST00000582591.1, ENST00000616852.1, MSTRG.11549.10  

### 3. Estimate parameters for sleuth reduced model
so_transcript <- sleuth_fit(so_transcript, ~1, "reduced")

### The models that have been fit can always be examined with the models() function:
# models(so_transcript)


########## Perform different types of analyses on the sleuth object

### 1. Perform differential analysis (testing) using the likelihood ratio test
so_transcript <- sleuth_lrt(so_transcript, "reduced", "full")

### 2. Perform differential analysis (testing) using the Wald test
so_transcript <- sleuth_wt(so_transcript, "conditionparental")


########## Examine results of test

q_val_threshold <- 0.025
sleuth_lrt_transcript <- sleuth_results(so_transcript, "reduced:full", "lrt", show_all = FALSE)
sleuth_lrt_transcript_significant <- filter(sleuth_lrt_transcript, qval <= q_val_threshold)
sleuth_wt_transcript <- sleuth_results(so_transcript, "conditionparental", show_all = TRUE)
sleuth_wt_transcript_significant <- filter(sleuth_wt_transcript, qval <= q_val_threshold)
# 7364

fold_threshold <- 2
# check how many genes were up- or down-regulated with a high/low 2-fold change
sleuth_wt_transcript_high_fold <- filter(sleuth_wt_transcript, abs(b) >= fold_threshold)
# 4639

########## Save into files

out_dir <- "../results/5_DiffExpr_3_DiffExpr_TransLevel"
write.csv(sleuth_lrt_transcript, file = paste0(out_dir, "/sleuth_lrt_transcript.csv"), row.names = FALSE)
write.csv(sleuth_lrt_transcript_significant, file = paste0(out_dir, "/sleuth_lrt_transcript_significant.csv"), row.names = FALSE)
write.csv(sleuth_wt_transcript, file = paste0(out_dir, "/sleuth_wt_transcript.csv"), row.names = FALSE)
write.csv(sleuth_wt_transcript_significant, file = paste0(out_dir, "/sleuth_wt_transcript_significant.csv"), row.names = FALSE)
sleuth_save(so_transcript, file = paste0(out_dir, "/so_transcript"))


########## Create plots

# Read in the sleuth object (only necessary if not done in same session and don't want to run everything again)
in_dir <- "../results/5_DiffExpr_3_DiffExpr_TransLevel"
so_transcript <- sleuth_load(paste0(in_dir, "/so_transcript"))

# Show sleuth live (tool for plots)
sleuth_live(so_transcript)


########## Read in the raw data (only needed to analyse if don't want to rerun the whole script)

# in_dir <- "../results/5_DiffExpr_3_DiffExpr_TransLevel"
# sleuth_lrt_transcript <- read.csv(file = paste0(in_dir, "/sleuth_lrt_transcript.csv"))
# sleuth_lrt_transcript_significant <- read.csv(file = paste0(in_dir, "/sleuth_lrt_transcript_significant.csv"))
# sleuth_wt_transcript <- read.csv(file = paste0(in_dir, "/sleuth_wt_transcript.csv"))
# sleuth_wt_transcript_significant <- read.csv(file = paste0(in_dir, "/sleuth_wt_transcript_significant.csv"))