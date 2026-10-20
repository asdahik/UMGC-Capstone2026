# Purpose: Import csv tables (Transcriptomics, Proteomics, Metabolomics, and Cultivar Features)
# Cleans and standardizes feature/phenotype IDs across all layers and export files for integration

# The following input files were downloaded as .csv from the Google Sheet data
# Transcriptomic_DESeq_PRJNA311702.csv : "Transcriptomics - DESeq PRJNA311702" tab from Google Sheet
# Proteomic_Abundance.csv : "Proteomic Abundance" tab from Google Sheet
# Cultivar_Metabolome.csv : "Cultivar Metabolome - Sprenger et al. 2016" tab from Google Sheet


# Load required libraries
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

wd <- getwd()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: gsea.r <PRJNA######>, <PXD######>")
PRJNA  <- args[1]
PXD <- args[2]
#META <- args [3] leaving meta argument out as we don't have any other metabolomic data at the moment.
trans_dir <- file.path(wd,"Transcriptomics","results",PRJNA)
prot_dir <- file.path(wd,"Proteomics", "results", PXD)
meta_dir <- file.path(wd, "Metabolomics", "cultivar_metabolome.txt") # this isn't a directory for right now

# Setting up directory and making an Integration folder if it isn't already there
wd <- getwd()
integration_dir <- file.path(wd, "Multiomics", "Integration")
if (!dir.exists(integration_dir)) dir.create(integration_dir, recursive = TRUE)


# Transcriptomics
# Reads in DESeq2 differential expression results from RNA seq scripts
# Creates a 4 column csv file with gene IDs, log2 fold change, adjusted p-value,
# and Omics type column for layer identification


# Transcriptomics data file from Google Sheet. Renamed when downloading
#rna_file <- "Transcriptomic_DESeq_PRJNA311702.csv"

# Read in the transcriptomics table exported from DESeq2 script
rna_raw <- read.csv(file.path(trans_dir, "deseq_results_annotated.csv"), stringsAsFactors = FALSE)


# Select and rename columns to standard 4-column format (ID, Log2FC, PValue, Omics)
rna_clean <- rna_raw %>%
  rename(
    ID     = SpudDBGeneID,
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
# prot_file <- "Proteomic_Abundance.csv"

# Read in the proteomics table
# Proteomics
prot_file <- file.path(prot_dir, paste0(PXD, "_Abundance.tsv"))

  # The header row has one extra leading name ("Protein FDR Confidence: Combined")
  # and data rows have one extra trailing empty field. Drop both so names line up.
  hdr <- strsplit(readLines(prot_file, n = 1), "\t", fixed = TRUE)[[1]]
  hdr <- trimws(hdr[-1])

  prot_raw <- read.delim(prot_file, header = FALSE, skip = 1, sep = "\t",
                        quote = "", comment.char = "", fill = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE)
  prot_raw <- prot_raw[, seq_along(hdr), drop = FALSE]   # drops trailing empty column
  colnames(prot_raw) <- hdr

  stress_col  <- "Abundance Ratio: (20% PEG, STRESS) / (CTRL, STRESS)"
  adapted_col <- "Abundance Ratio: (20% PEG, Adapted) / (CTRL, Adapted)"
  pval_col    <- "Abundance Ratio Adj. P-Value: (20% PEG, STRESS) / (CTRL, STRESS)"

  missing_cols <- setdiff(c("Accession", stress_col, adapted_col, pval_col), colnames(prot_raw))
  if (length(missing_cols) > 0) stop("Missing columns:\n", paste(missing_cols, collapse = "\n"))

  prot_clean <- tibble(
    ID            = trimws(prot_raw[["Accession"]]),
    Stress_Ratio  = suppressWarnings(as.numeric(prot_raw[[stress_col]])),
    Adapted_Ratio = suppressWarnings(as.numeric(prot_raw[[adapted_col]])),
    PValue        = suppressWarnings(as.numeric(prot_raw[[pval_col]]))   # already 0-1
  ) %>%
    mutate(
      Log2FC = log2((Stress_Ratio + 1e-5) / (Adapted_Ratio + 1e-5)),
      Omics  = "Protein"
    ) %>%
    filter(!is.na(ID), ID != "", !is.na(Log2FC), !is.infinite(Log2FC)) %>%
    select(ID, Log2FC, PValue, Omics)

  write.csv(prot_clean, file.path(integration_dir, "proteomics_clean.csv"), row.names = FALSE)


# Metabolomics
  
# Metabolomics data file from Google Sheet. Renamed when downloading
#metab_file <- "Cultivar_Metabolome.csv"

  # Read in the metabolomics table. Skip the first 12  description rows 
  # so row 13 acts as the true header
# Tab-delimited; rows 1-12 are description rows, row 13 is the true header
metab_raw <- read.delim(meta_dir, skip = 12, header = TRUE, sep = "\t",
                        quote = "", comment.char = "", check.names = FALSE,
                        stringsAsFactors = FALSE, na.strings = c("NA", ""))

  stopifnot("Treatment" %in% colnames(metab_raw))   # fails with a clear message if misparsed
  message("Metabolomics dims: ", nrow(metab_raw), " x ", ncol(metab_raw))
  print(table(tolower(trimws(metab_raw$Treatment)), useNA = "ifany"))

  # Pull treatment labels and make lowercase
  treatment_col <- tolower(trimws(metab_raw[["Treatment"]]))
  # Isolate just the metabolite abundance columns (starting at column 11)
  metab_matrix  <- metab_raw[, 11:ncol(metab_raw), drop = FALSE]
  
  
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

hdr_lines <- readLines(meta_dir, n = 12)
row_vals  <- function(i) strsplit(hdr_lines[i], "\t", fixed = TRUE)[[1]]
idx       <- 11:(10 + ncol(metab_matrix))   # descriptor values start in column 11, like the data

clean_na <- function(x) { x <- gsub('^"|"$', "", trimws(x)); x[x == "NA" | x == ""] <- NA; x }

metab_map <- data.frame(
  ID     = colnames(metab_matrix),
  Name   = clean_na(row_vals(4)[idx]),
  Class  = clean_na(row_vals(3)[idx]),
  KEGG   = clean_na(row_vals(7)[idx]),
  GMD_ID = clean_na(row_vals(11)[idx]),
  stringsAsFactors = FALSE
)

# Where the name row just repeats the code, the analyte is unidentified
metab_map$Annotated <- !is.na(metab_map$Name) & metab_map$Name != metab_map$ID

# Sanity check: A204001-101 should show its name and KEGG ID here
print(metab_map[metab_map$ID == "A204001-101", ])

write.csv(metab_map, file.path(integration_dir, "metabolomics_annotation.csv"), row.names = FALSE)