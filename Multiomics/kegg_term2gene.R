# Purpose: Build KEGG TERM2GENE / TERM2NAME tables from an eggNOG-mapper annotation file
# Usage:   Rscript build_kegg_term2gene.R <emapper.annotations> <out_prefix> [strip_isoform: TRUE/FALSE]
#   e.g.   Rscript build_kegg_term2gene.R v61_proteome.emapper.annotations Transcriptomics/ref/term2gene_KEGG TRUE
#
# Writes (no header, tab-separated, same layout as term2gene_GO_BP.tsv):
#   <out_prefix>.tsv          term (ko#####)  <tab>  gene
#   <out_prefix_name>.tsv     term (ko#####)  <tab>  pathway name
#
# For the RNA layer, queries are transcript/protein IDs like Soltu.DM.06G030080.1,
# so strip_isoform = TRUE collapses them to the gene ID Soltu.DM.06G030080.
# For the proteomics layer (UniProt accessions), use strip_isoform = FALSE.

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: build_kegg_term2gene.R <emapper.annotations> <out_prefix> [strip_isoform]")
ann_file       <- args[1]
out_prefix     <- args[2]
strip_isoform  <- if (length(args) >= 3) as.logical(args[3]) else TRUE

# eggNOG-mapper files have '##' comment lines; the header line starts with '#query'
lines <- readLines(ann_file, warn = FALSE)
lines <- lines[!grepl("^##", lines)]
ann   <- read.delim(text = lines, sep = "\t", quote = "", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE)
colnames(ann)[1] <- "query"
stopifnot("KEGG_Pathway" %in% colnames(ann))

message("Annotated queries: ", nrow(ann),
        " | with KEGG pathway: ", sum(ann$KEGG_Pathway != "-"))

# KEGG_Pathway format differs by eggNOG-mapper version: newer releases give bare
# 5-digit codes ("00220,00250,04136"); older ones give "ko00010,map00010".
# Normalise everything to ko##### and keep each pathway once per gene.
t2g <- ann %>%
  transmute(gene = if (strip_isoform) sub("\\.\\d+$", "", query) else query,
            term = KEGG_Pathway) %>%
  filter(term != "-", !is.na(term)) %>%
  separate_rows(term, sep = ",") %>%
  mutate(term = sub("^(ko|map)", "", trimws(term))) %>%
  filter(grepl("^\\d{5}$", term)) %>%
  mutate(term = paste0("ko", term)) %>%
  filter(!grepl("^ko01\\d{3}$", term)) %>%    # drop global/overview maps (01100, 01110, 01200, ...)
  distinct(term, gene)

# Pathway names from KEGG (needs internet). Column layout can differ by
# clusterProfiler version, so check head(nm) if the join below looks empty.
kegg <- download_KEGG("ko")
nm   <- kegg$KEGGPATHID2NAME
colnames(nm) <- c("term", "name")
nm$term <- ifelse(grepl("^ko", nm$term), nm$term, paste0("ko", sub("^[a-z]+", "", nm$term)))
t2n <- nm %>% filter(term %in% t2g$term) %>% distinct()

# KO-based assignment also hits pathways that don't exist in plants (e.g. human
# disease maps). Optionally keep only pathways KEGG lists for potato (organism "sot").
plant_ids <- tryCatch({
  sot <- download_KEGG("sot")$KEGGPATHID2EXTID
  unique(paste0("ko", sub("^[a-z]+", "", sot[[1]])))
}, error = function(e) { message("Could not fetch sot pathway list; keeping all pathways"); NULL })
if (!is.null(plant_ids)) {
  t2g <- t2g %>% filter(term %in% plant_ids)
  t2n <- t2n %>% filter(term %in% plant_ids)
}

message("Pathways: ", n_distinct(t2g$term), " | genes: ", n_distinct(t2g$gene),
        " | pairs: ", nrow(t2g), " | named pathways: ", nrow(t2n))
print(head(t2g, 3))
print(head(t2n, 3))

write.table(t2g[, c("term", "gene")], paste0(out_prefix, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(t2n, paste0(out_prefix, "_names.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)