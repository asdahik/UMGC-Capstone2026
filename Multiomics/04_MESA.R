# Purpose: Metabolite set enrichment on KEGG pathways (GSEA-style + ORA), then an
#          optional pathway-level table lining up RNA, protein and metabolite results.
# Usage:   Rscript metabolomics_gsea.R [PRJNA######] [PXD######]
#          The two arguments are only needed for the cross-layer table at the end.
#
# Needs:   Metabolomics/cultivar_metabolome.txt                (tab-delimited, 12 description rows, row 13 = header)
#          Transcriptomics/ref/term2gene_KEGG_names.tsv        (the 141 potato pathways, ko##### + name)
#          Internet access the first time (KEGG compound->pathway links; cached afterwards)
#
# Pathway IDs are ko##### everywhere so the layers can be merged by ID.
 
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(KEGGREST)
})
 
args  <- commandArgs(trailingOnly = TRUE)
PRJNA <- if (length(args) >= 1) args[1] else NA_character_
PXD   <- if (length(args) >= 2) args[2] else NA_character_
 
wd        <- getwd()
meta_file <- file.path(wd, "Metabolomics", "cultivar_metabolome.txt")
ref_dir   <- file.path(wd, "Transcriptomics", "ref")
mo_dir    <- file.path(wd, "Multiomics", "results")
out_dir   <- file.path(wd, "Multiomics", "Integration", "GSEA_metabolomics")
dir.create(mo_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
 
# ---- settings you may want to change -----------------------------------
min_set      <- 3      # smallest pathway (number of measured metabolites) to test
max_set      <- 500
ora_p_cutoff <- 0.05   # for ORA, applied to BH-adjusted p-values
min_ora_n    <- 3      # skip ORA lists with fewer metabolites than this
# Optional subsetting before the t-test. Leave NULL to use every non-control sample.
stress_levels <- NULL  # e.g. c("drought")  (lower-case Treatment labels)
trial_subset  <- NULL  # e.g. c("A1", "A2")  (values of "Trial ID")
 
# ---- 1. Read the metabolome ---------------------------------------------
# Rows 1-12 describe each analyte (row 3 Class, 4 Name, 7 KEGG-ID, 11 GMD-ID);
# row 13 is the true header. Analyte columns start at column 11.
n_col <- max(count.fields(meta_file, sep = "\t", quote = "\"", comment.char = "")[1:13],
             na.rm = TRUE)
desc <- read.delim(meta_file, header = FALSE, nrows = 12, sep = "\t", quote = "\"",
                   comment.char = "", fill = TRUE, colClasses = "character",
                   col.names = paste0("V", seq_len(n_col)), na.strings = character(0))
 
metab_raw <- read.delim(meta_file, skip = 12, header = TRUE, sep = "\t", quote = "",
                        comment.char = "", check.names = FALSE,
                        na.strings = c("NA", ""), stringsAsFactors = FALSE)
stopifnot("Treatment" %in% colnames(metab_raw))
 
met_idx  <- 11:ncol(metab_raw)
met_cols <- colnames(metab_raw)[met_idx]
clean_na <- function(x) { x <- trimws(x); x[x %in% c("NA", "")] <- NA; x }
pull_row <- function(r) clean_na(unlist(desc[r, met_idx], use.names = FALSE))
 
met_map <- tibble(ID = met_cols, Name = pull_row(4), Class = pull_row(3),
                  KEGG = pull_row(7))
 
# Sanity check: where the Name row just repeats an analyte code, it must equal the header
code_like <- !is.na(met_map$Name) & grepl("^A\\d{6}-101$", met_map$Name)
if (any(code_like) && !all(met_map$Name[code_like] == met_map$ID[code_like]))
  warning("Descriptor rows do not line up with the data header; check the 12-row layout.")
 
message("Analytes: ", nrow(met_map), " | with KEGG compound ID: ", sum(!is.na(met_map$KEGG)))
print(table(Treatment = tolower(trimws(metab_raw$Treatment)), useNA = "ifany"))
 
# ---- 2. Per-metabolite statistics (stress vs control) --------------------
trt  <- tolower(trimws(metab_raw$Treatment))
keep <- rep(TRUE, nrow(metab_raw))
if (!is.null(trial_subset)) keep <- keep & metab_raw[["Trial ID"]] %in% trial_subset
is_ctrl   <- keep & !is.na(trt) & trt == "control"
is_stress <- keep & !is.na(trt) & trt != "" & trt != "control"
if (!is.null(stress_levels)) is_stress <- is_stress & trt %in% tolower(stress_levels)
message("Control samples: ", sum(is_ctrl), " | stress samples: ", sum(is_stress))
 
stats <- lapply(met_cols, function(cn) {
  v  <- suppressWarnings(as.numeric(metab_raw[[cn]]))
  s  <- v[is_stress & !is.na(v)]
  c0 <- v[is_ctrl   & !is.na(v)]
  if (length(s) < 2 || length(c0) < 2) return(NULL)
  tt <- tryCatch(t.test(s, c0), error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  # values are log10 intensities, so the mean difference * log2(10) is a log2 fold change
  tibble(ID = cn, Log2FC = (mean(s) - mean(c0)) * log2(10),
         tstat = unname(tt$statistic), PValue = tt$p.value)
}) %>% bind_rows()
 
stats <- stats %>% mutate(padj = p.adjust(PValue, method = "BH")) %>% left_join(met_map, by = "ID")
write.csv(stats, file.path(out_dir, "metabolite_stats.csv"), row.names = FALSE)
message("Metabolites tested: ", nrow(stats), " | BH-adjusted p < 0.05: ", sum(stats$padj < 0.05, na.rm = TRUE))
 
# Keep identified metabolites; if two analytes share a KEGG ID keep the stronger one
kstats <- stats %>%
  filter(!is.na(KEGG), is.finite(tstat)) %>%
  group_by(KEGG) %>% slice_max(abs(tstat), n = 1, with_ties = FALSE) %>% ungroup()
message("Metabolites usable for pathway analysis: ", nrow(kstats))
 
# ---- 3. Compound -> pathway sets (potato pathways only, ko##### IDs) ----
k2n_raw <- read.delim(file.path(ref_dir, "term2gene_KEGG_names.tsv"), header = FALSE,
                      col.names = c("term", "name"), quote = "", stringsAsFactors = FALSE)
k2n <- k2n_raw %>% mutate(name = paste0(name, " (", term, ")"))
 
cache <- file.path(ref_dir, "compound2pathway_KEGG.tsv")
if (file.exists(cache)) {
  c2p <- read.delim(cache, header = FALSE, col.names = c("term", "compound"),
                    stringsAsFactors = FALSE)
} else {
  lk  <- keggLink("pathway", "compound")
  c2p <- tibble(compound = sub("^cpd:", "", names(lk)), path = sub("^path:", "", unname(lk))) %>%
    filter(grepl("^(map|ko)\\d{5}$", path)) %>%
    mutate(term = paste0("ko", sub("^[a-z]+", "", path))) %>%
    distinct(term, compound)
  write.table(c2p[, c("term", "compound")], cache, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}
c2p <- c2p %>% filter(term %in% k2n$term, compound %in% kstats$KEGG) %>%
  rename(gene = compound)
sizes <- c2p %>% count(term, name = "n")
message("Pathways with >= ", min_set, " measured metabolites: ", sum(sizes$n >= min_set),
        " of ", nrow(sizes), " with any")
 
# ---- 4. GSEA-style enrichment (ranked by t-statistic) --------------------
ranks <- setNames(kstats$tstat, kstats$KEGG)
ranks <- sort(ranks, decreasing = TRUE)
 
gsea_m <- tryCatch(
  GSEA(ranks, TERM2GENE = c2p, TERM2NAME = k2n, minGSSize = min_set, maxGSSize = max_set,
       pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE),
  error = function(e) { message("Metabolite GSEA failed: ", conditionMessage(e)); NULL })
res_m <- if (is.null(gsea_m)) data.frame() else as.data.frame(gsea_m)
 
if (nrow(res_m) > 0) {
  cat("\nMetabolite GSEA: ", nrow(res_m), " pathways tested, ",
      sum(res_m$p.adjust < 0.05), " at BH < 0.05, ",
      sum(res_m$pvalue < 0.05), " at raw p < 0.05\n", sep = "")
  print(head(res_m[order(res_m$pvalue), c("Description", "setSize", "NES", "pvalue", "p.adjust")], 15))
  write.csv(res_m, file.path(out_dir, "metabolite_GSEA_KEGG.csv"), row.names = FALSE)
  write.csv(res_m, file.path(mo_dir, "metabolite_gsea_KEGG.csv"), row.names = FALSE)
} else {
  message("No pathways passed the size filter for GSEA.")
}
 
# ---- 5. ORA cross-check on significant metabolites -----------------------
universe <- kstats$KEGG
sig <- kstats %>% filter(padj < ora_p_cutoff)
lists <- list(all_changed = sig$KEGG,
              up          = sig$KEGG[sig$tstat > 0],
              down        = sig$KEGG[sig$tstat < 0])
 
ora_all <- list()
for (nm_ in names(lists)) {
  g <- lists[[nm_]]
  if (length(g) < min_ora_n) {
    message("ORA '", nm_, "': only ", length(g), " significant metabolites; skipped.")
    next
  }
  r <- tryCatch(
    enricher(g, universe = universe, TERM2GENE = c2p, TERM2NAME = k2n,
             minGSSize = min_set, maxGSSize = max_set, pvalueCutoff = 1, qvalueCutoff = 1),
    error = function(e) NULL)
  if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
    ora_all[[nm_]] <- as.data.frame(r) %>% mutate(set = nm_)
    cat("\nORA '", nm_, "' (", length(g), " metabolites):\n", sep = "")
    print(head(ora_all[[nm_]][, c("Description", "Count", "GeneRatio", "pvalue", "p.adjust")], 10))
  }
}
if (length(ora_all) > 0)
  write.csv(bind_rows(ora_all), file.path(out_dir, "metabolite_ORA_KEGG.csv"), row.names = FALSE)
 
# Pathways flagged by both approaches (raw p < 0.05 in each; sets are tiny, so adjusted
# p-values are rarely informative here)
if (nrow(res_m) > 0 && length(ora_all) > 0) {
  both <- intersect(res_m$ID[res_m$pvalue < 0.05],
                    bind_rows(ora_all)$ID[bind_rows(ora_all)$pvalue < 0.05])
  cat("\nPathways with raw p < 0.05 in both GSEA and ORA: ", length(both), "\n", sep = "")
  if (length(both) > 0) print(k2n[k2n$term %in% both, ])
}
 
# ---- 6. Cross-layer pathway table (RNA + protein + metabolite) -----------
norm_ko <- function(x) paste0("ko", sub("^[A-Za-z]+", "", x))
read_layer <- function(path, tag) {
  if (is.na(path) || !file.exists(path)) {
    message("Layer '", tag, "' not found, skipped: ", path); return(NULL)
  }
  d <- read.csv(path, stringsAsFactors = FALSE)
  tibble(ID = norm_ko(d$ID), NES = d$NES, padj = d$p.adjust, pval = d$pvalue) %>%
    rename_with(~ paste0(.x, "_", tag), -ID)
}
layers <- list(
  read_layer(if (!is.na(PRJNA)) file.path(mo_dir, paste0(PRJNA, "_gsea_KEGG.csv")) else NA, "RNA"),
  read_layer(file.path(wd, "Multiomics", "results",
                       "proteomics_GSEA_KEGG.csv"), "Protein"),
  read_layer(file.path(mo_dir, "metabolite_gsea_KEGG.csv"), "Metabolite")
)
layers <- Filter(Negate(is.null), layers)
 
if (length(layers) >= 2) {
  integ <- Reduce(function(a, b) full_join(a, b, by = "ID"), layers) %>%
    left_join(k2n_raw %>% rename(ID = term, Pathway = name), by = "ID")
  nes_cols  <- grep("^NES_",  colnames(integ), value = TRUE)
  padj_cols <- grep("^padj_", colnames(integ), value = TRUE)
  integ$n_layers      <- rowSums(!is.na(integ[, nes_cols, drop = FALSE]))
  integ$n_sig         <- rowSums(integ[, padj_cols, drop = FALSE] < 0.05, na.rm = TRUE)
  integ$same_sign_all <- apply(integ[, nes_cols, drop = FALSE], 1, function(x) {
    x <- x[!is.na(x)]; length(x) >= 2 && (all(x > 0) || all(x < 0))
  })
  integ <- integ %>% arrange(desc(n_sig), desc(same_sign_all), desc(n_layers))
  write.csv(integ, file.path(mo_dir, "KEGG_pathway_integration.csv"), row.names = FALSE)
  cat("\nCross-layer table written: ", nrow(integ), " pathways, ",
      sum(integ$n_layers >= 2 & integ$same_sign_all, na.rm = TRUE),
      " with the same direction in 2+ layers\n", sep = "")
  print(head(integ %>% select(Pathway, starts_with("NES_"), n_sig), 15))
} else {
  message("Fewer than two layers available; cross-layer table not written.")
}
 
message("\nDone. Outputs in: ", out_dir)
# test correlation between found values
ok <- complete.cases(integ$NES_RNA, integ$NES_Protein)
cor.test(integ$NES_RNA[ok], integ$NES_Protein[ok], method = "spearman")
plot(integ$NES_RNA[ok], integ$NES_Protein[ok])   # or ggplot with labels for the significant one