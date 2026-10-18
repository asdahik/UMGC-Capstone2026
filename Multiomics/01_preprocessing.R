# Purpose: Import csv tables (Transcriptomics, Proteomics, Metabolomics, and Cultivar Features)
# Cleans and standardizes feature/phenotype IDs across all layers and export files for integration

# The following input files were downloaded as .csv from the Google Sheet data
# Transcriptomic_DESeq_PRJNA311702.csv : "Transcriptomics - DESeq PRJNA311702" tab from Google Sheet
# Proteomic_Abundance.csv : "Proteomic Abundance" tab from Google Sheet
# Cultivar_Metabolome.csv : "Cultivar Metabolome - Sprenger et al. 2016" tab from Google Sheet


# Load required libraries
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# Setting up directory and making an Integration folder if it isn't already there
wd <- getwd()
integration_dir <- file.path(wd, "Integration")
if (!dir.exists(integration_dir)) dir.create(integration_dir, recursive = TRUE)


# Transcriptomics
# Reads in DESeq2 differential expression results from RNA seq scripts
# Creates a 4 column csv file with gene IDs, log2 fold change, adjusted p-value,
# and Omics type column for layer identification


# Transcriptomics data file from Google Sheet. Renamed when downloading
rna_file <- "Transcriptomic_DESeq_PRJNA311702.csv"


# Read in the transcriptomics table exported from DESeq2 script
rna_raw <- read.csv(rna_file, stringsAsFactors = FALSE)


# Select and rename columns to standard 4-column format (ID, Log2FC, PValue, Omics)
rna_clean <- rna_raw %>%
  rename(
    ID     = SpudDB.Gene.ID,
    Log2FC = log2FoldChange,
    PValue = padj
  ) %>%
  filter(!is.na(ID), ID != "", !is.na(Log2FC)) %>%
  mutate(Omics = "RNA") %>%
  select(ID, Log2FC, PValue, Omics)

# Save the cleaned table into integration folder
write.csv(rna_clean, file.path(integration_dir, "transcriptomics_clean.csv"), row.names = FALSE)



# Proteomics 

# Proteomics data file from Google Sheet 
prot_file <- "Proteomic_Abundance.csv"

# Read in the proteomics table
prot_raw <- read.csv(prot_file, stringsAsFactors = FALSE, check.names = FALSE)

# Pull out protein IDs, calculate log2 fold changes between stress and adapted
# (added 1e-5 to avoid division by zero errors),
# and scale the percentage based p-values down to 0-1 range
  prot_clean <- prot_raw %>%
    mutate(
      ID            = prot_raw[[1]], 
      Stress_Ratio  = suppressWarnings(as.numeric(`Abundance Ratio: (20% PEG, STRESS) / (CTRL, STRESS)`)),
      Adapted_Ratio = suppressWarnings(as.numeric(`Abundance Ratio: (20% PEG, Adapted) / (CTRL, Adapted)`)),
      Log2FC        = log2((Stress_Ratio + 1e-5) / (Adapted_Ratio + 1e-5)),
      # Divide by 100 to convert the percentage-scaled column into a true 0 to 1 p-value
      PValue        = suppressWarnings(as.numeric(`Abundance Ratio Adj. P-Value: (20% PEG, STRESS) / (CTRL, STRESS)`)) / 100,
      Omics         = "Protein"
    ) %>%
    filter(!is.na(ID), ID != "", !is.na(Log2FC), !is.infinite(Log2FC)) %>%
    select(ID, Log2FC, PValue, Omics)

  # Save the cleaned table into integration folder
  write.csv(prot_clean, file.path(integration_dir, "proteomics_clean.csv"), row.names = FALSE)


# Metabolomics
  
# Metabolomics data file from Google Sheet. Renamed when downloading
metab_file <- "Cultivar_Metabolome.csv"

  # Read in the metabolomics table. Skip the first 12  description rows 
  # so row 13 acts as the true header
  metab_raw <- read.csv(metab_file, check.names = FALSE, skip = 12, stringsAsFactors = FALSE)
  
  # Pull treatment labels and make lowercase
  treatment_col <- tolower(metab_raw[["Treatment"]])
  
  # Isolate just the metabolite abundance columns (starting at column 11)
  metab_matrix  <- metab_raw[, 11:ncol(metab_raw)] 
  
  # Initialize vectors to store results
  metabolite_ids <- c()
  log2fc_vals    <- c()
  p_vals         <- c()
  
  for (col_name in colnames(metab_matrix)) {
    # Convert column values to numeric
    vals <- suppressWarnings(as.numeric(as.character(metab_matrix[[col_name]])))
    
    # Split by treatment
    ctrl_vals   <- vals[treatment_col == "control"]
    stress_vals <- vals[treatment_col != "control"]
    
    # Remove NAs for calculation
    ctrl_clean   <- ctrl_vals[!is.na(ctrl_vals)]
    stress_clean <- stress_vals[!is.na(stress_vals)]
    
    # Ensure we have enough replicates to calculate statistics
    if (length(ctrl_clean) > 1 && length(stress_clean) > 1) {
      ctrl_mean   <- mean(ctrl_clean)
      stress_mean <- mean(stress_clean)
      
      # Calculate Log2 Fold Change (Note: source values are log10 transformed in this table)
      # Log2FC = (Mean_Stress - Mean_Control) * log2(10) to properly scale to base 2
      fc_log2 <- (stress_mean - ctrl_mean) * log2(10)
      
      # Run t-test to get a true statistical p-value
      t_test <- t.test(stress_clean, ctrl_clean)
      
      # Add results to vectors
      metabolite_ids <- c(metabolite_ids, col_name)
      log2fc_vals    <- c(log2fc_vals, fc_log2)
      p_vals         <- c(p_vals, t_test$p.value)
    }
  }
  
  # Build clean metabolomics dataframe
  metab_clean <- data.frame(
    ID     = metabolite_ids,
    Log2FC = log2fc_vals,
    PValue = p_vals,
    Omics  = "Metabolite",
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(ID), ID != "", !is.na(Log2FC), !is.infinite(Log2FC))
 
  #  Save the cleaned table into integration folder
  write.csv(metab_clean, file.path(integration_dir, "metabolomics_clean.csv"), row.names = FALSE)

