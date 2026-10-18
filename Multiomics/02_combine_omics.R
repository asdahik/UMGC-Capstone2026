# Purpose: Combines transcriptomics, proteomics, and metabolomics data.
# Adds regulation labels without dropping any features.
# Feature selection will occur after predictive modeling

library(dplyr)

# Set up integration directory path
integration_dir <- file.path(getwd(), "Integration")

# Read in the clean individual layers exported from 01_preprocessing.R
rna_clean <- read.csv(file.path(integration_dir, "transcriptomics_clean.csv"), stringsAsFactors = FALSE)
prot_clean <- read.csv(file.path(integration_dir, "proteomics_clean.csv"), stringsAsFactors = FALSE)
metab_clean <- read.csv(file.path(integration_dir, "metabolomics_clean.csv"), stringsAsFactors = FALSE)

# Combine all three clean tables vertically
multi_omics <- bind_rows(rna_clean, prot_clean, metab_clean)

# Save the multi omics file in integration folder
write.csv(multi_omics, file.path(integration_dir, "multi_omics_clean.csv"), row.names = FALSE)
