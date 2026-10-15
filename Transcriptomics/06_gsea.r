library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(GO.db)

wd <- getwd()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: gsea.r <PRJNA######>")
PRJNA  <- args[1]
res_dir <- file.path(wd,"Transcriptomics","results",PRJNA)
counts_file  <- file.path(res_dir, paste0(PRJNA, "_counts.txt"))   # paste0 instead of + (+ cannot join strings in R)
meta_file <- read.csv(file.path(res_dir,"SraRunTable.csv"))
tab <- read.csv(file.path(res_dir, "deseq_results_annotated.csv"), row.names = 1)

fc <- read.delim(counts_file, comment.char ="#", check.names = FALSE)

count_mat <- fc[, 7:ncol(fc)]  # Columns 7+ are the sample count_mat
rownames(count_mat) <- fc$Geneid
colnames(count_mat) <- sub("\\.sorted\\.bam$", "", basename(colnames(count_mat)))

# meta data for sampling metadata
meta <- meta_file[match(colnames(count_mat), meta_file$Run), ] #column run stores SRA
stopifnot(all(meta$Run == colnames(count_mat)))
rownames(meta) <- meta$Run

meta$treatment <- factor(sub(" ", "_", meta$treatment), levels = c("control", "drought_stress"))
meta$Cultivar <- factor(meta$Cultivar)
meta$experiment <- factor(meta$experiment)
# From the GEO growth protocol: F1/F3/F4 = field, G1/G2/G3 = greenhouse
meta$environment <- factor(ifelse(grepl("^F", as.character(meta$experiment)),
                                  "field", "greenhouse"))
print(table(meta$experiment, meta$environment))

# pooled ranking: Wald statistic from your pooled drought-vs-control results
ranks <- tab$log2FoldChange / tab$lfcSE          # Wald statistic
names(ranks) <- rownames(tab)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

# experiment x treatment interaction test (informational only; its statistic
# is NOT used for GSEA because the LRT statistic has no up/down direction)
dds_exp <- DESeqDataSetFromMatrix(count_mat, meta,
             design = ~ experiment + Cultivar + treatment + experiment:treatment)
keep <- rowSums(counts(dds_exp)) >= 10
dds_exp <- dds_exp[keep, ]
dds_exp <- DESeq(dds_exp, test = "LRT",
                 reduced = ~ experiment + Cultivar + treatment)
res_exp <- results(dds_exp)
cat(sum(res_exp$padj < 0.05, na.rm = TRUE), "genes with a significant experiment x treatment interaction\n")

# create tables for
t2g <- read.delim(file.path(wd,"Transcriptomics","ref","term2gene_GO_BP.tsv"), header = FALSE,
                  col.names = c("term", "gene"), stringsAsFactors = FALSE)

# drop obsolete / unrecognized GO terms (terms without a name in GO.db)
nm0 <- AnnotationDbi::select(GO.db, keys = unique(t2g$term),
                             columns = "TERM", keytype = "GOID")
t2g <- t2g[t2g$term %in% nm0$GOID[!is.na(nm0$TERM)], ]

# filter input from term to gene table
t2g <- t2g[t2g$gene %in% names(ranks), ]

# This step is for getting ancestors for each gene
# ancestors of each BP term (as.list gives a named list)
anc <- as.list(GOBPANCESTOR)
terms <- unique(t2g$term)
anc <- anc[names(anc) %in% terms]

anc_tab <- data.frame(child  = rep(names(anc), lengths(anc)),
                      parent = unlist(anc, use.names = FALSE),
                      stringsAsFactors = FALSE)
anc_tab <- anc_tab[anc_tab$parent != "all", ]               # drop the root placeholder
anc_tab <- rbind(anc_tab, data.frame(child = terms, parent = terms))   # keep each term itself

# expand term to gene with the ancestor table
t2g_exp <- unique(merge(t2g, anc_tab, by.x = "term", by.y = "child")[, c("parent", "gene")])
colnames(t2g_exp) <- c("term", "gene")

# readable term names
nm <- AnnotationDbi::select(GO.db, keys = unique(t2g_exp$term),
                            columns = "TERM", keytype = "GOID")
t2n <- data.frame(term = nm$GOID, name = nm$TERM, stringsAsFactors = FALSE)

cat(length(ranks), "ranked genes\n")

cat(nrow(t2g), "pairs before,", nrow(t2g_exp), "after expansion\n")
write.table(t2g_exp, "term2gene_GO_BP_expanded.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(t2n, "term2name_GO_BP.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# change term to name varible to the new written table
t2n <- read.delim(file.path(wd,"Transcriptomics","ref","term2name_GO_BP.tsv"), header = FALSE,
                  col.names = c("term", "name"), quote = "", stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# Pooled GSEA (your block, with the result assigned and t2g_exp used)
# ------------------------------------------------------------------
g <- t2g_exp[t2g_exp$gene %in% names(ranks), ]
gsea <- GSEA(ranks, TERM2GENE = g, TERM2NAME = t2n,
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
             eps = 0, seed = TRUE, verbose = FALSE)
res_g <- as.data.frame(gsea)

cat(sum(res_g$p.adjust < 0.05), "enriched terms at p.adjust < 0.05\n")
print(head(res_g[order(res_g$p.adjust), c("Description", "setSize", "NES", "p.adjust")], 20))

write.csv(res_g, file.path(res_dir, "gsea_GO_BP.csv")), row.names = FALSE)

p <- dotplot(gsea, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign)
ggsave(file.path(res_dir, paste0(PRJNA, "gsea_GO_BP_dotplot.pdf")), p, width = 10, height = 8)

# use only gene pathways that have a biological process
res_g <- read.csv(file.path(res_dir,"gsea_GO_BP.csv"))
res_g <- res_g[order(-res_g$NES), ]
cols <- c("Description", "setSize", "NES", "p.adjust")

cat(sum(res_g$NES > 0 & res_g$p.adjust < 0.05), "terms up,",
    sum(res_g$NES < 0 & res_g$p.adjust < 0.05), "terms down\n")
cat("Strongest Up\n")
print(head(res_g[, cols], 15))            # strongest up
cat("Strongest Down\n\n")
print(tail(res_g[, cols], 15))            # strongest down

# ------------------------------------------------------------------
# Field and greenhouse: separate DESeq2 fits
# ------------------------------------------------------------------
fit_deseq <- function(idx) {
  m <- droplevels(meta[idx, ])
  d <- DESeqDataSetFromMatrix(count_mat[, idx], m,
                              design = ~ experiment + Cultivar + treatment)
  d <- d[rowSums(counts(d)) >= 10, ]
  d <- DESeq(d)
  results(d, contrast = c("treatment", "drought_stress", "control"), alpha = 0.05)
}

res_field <- fit_deseq(meta$environment == "field")
res_gh    <- fit_deseq(meta$environment == "greenhouse")

saveRDS(list(field = res_field, greenhouse = res_gh),
        file.path(res_dir, "res_by_environment.rds"))

cat("\n--- field ---\n");      summary(res_field)
cat("\n--- greenhouse ---\n"); summary(res_gh)

# gene-level agreement between the environments
common <- intersect(rownames(res_field), rownames(res_gh))
rf <- res_field[common, ]
rg <- res_gh[common, ]

cat("\nSpearman correlation of Wald statistics:",
    round(cor(rf$stat, rg$stat, use = "complete.obs", method = "spearman"), 3), "\n")
print(table(field = rf$padj < 0.05, greenhouse = rg$padj < 0.05, useNA = "ifany"))

ok    <- !is.na(rf$padj) & !is.na(rg$padj)
both  <- ok & rf$padj < 0.05 & rg$padj < 0.05
agree <- sign(rf$log2FoldChange) == sign(rg$log2FoldChange)
cat(sum(both), "genes significant in both;",
    round(100 * mean(agree[both]), 1), "% agree in direction\n")

shared_up   <- common[both & agree & rf$log2FoldChange > 0]
shared_down <- common[both & agree & rf$log2FoldChange < 0]
discordant  <- common[both & !agree]
universe    <- common[ok]
cat(length(shared_up), "shared up,", length(shared_down), "shared down,",
    length(discordant), "discordant\n")

for (nm_ in c("shared_up", "shared_down", "discordant")) {
  write.table(get(nm_), file.path(res_dir, paste0(nm_, "_genes.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ------------------------------------------------------------------
# GSEA per environment, then the combined comparison
# ------------------------------------------------------------------
run_gsea <- function(res) {
  rk <- res$stat
  names(rk) <- rownames(res)
  rk <- sort(rk[is.finite(rk)], decreasing = TRUE)
  gs <- t2g_exp[t2g_exp$gene %in% names(rk), ]
  as.data.frame(GSEA(rk, TERM2GENE = gs, TERM2NAME = t2n,
                     minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
                     eps = 0, seed = TRUE, verbose = FALSE))
}

gf <- run_gsea(res_field)
gg <- run_gsea(res_gh)
gp <- read.csv(file.path(res_dir, "gsea_GO_BP.csv"))      # pooled, written above

pick <- function(x, suffix) {
  x <- x[, c("ID", "NES", "p.adjust")]
  colnames(x)[2:3] <- paste0(c("NES_", "padj_"), suffix)
  x
}
cmp <- merge(gf[, c("ID", "Description", "setSize")], pick(gf, "field"), by = "ID")
cmp <- merge(cmp, pick(gg, "gh"),     by = "ID")
cmp <- merge(cmp, pick(gp, "pooled"), by = "ID")

sf   <- cmp$padj_field < 0.05
sg   <- cmp$padj_gh    < 0.05
same <- sign(cmp$NES_field) == sign(cmp$NES_gh)
cmp$class <- ifelse(sf & sg & same, "both, same direction",
             ifelse(sf & sg,        "both, opposite direction",
             ifelse(sf,             "field only",
             ifelse(sg,             "greenhouse only", "neither"))))

cat("\nNES correlation field vs greenhouse:",
    round(cor(cmp$NES_field, cmp$NES_gh), 3), "\n")
print(table(cmp$class))

write.csv(cmp, file.path(res_dir, paste0(PRJNA, "_gsea.csv")), row.names = FALSE)

p2 <- ggplot(cmp, aes(NES_field, NES_gh, color = class)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  theme_bw()
ggsave(file.path(res_dir, paste0(PRJNA, "_gsea_comparison.pdf")), p2, width = 8, height = 6)

shared <- cmp[cmp$class == "both, same direction", ]
shared$min_abs_NES <- pmin(abs(shared$NES_field), abs(shared$NES_gh))
shared <- shared[order(-shared$min_abs_NES), ]
write.csv(shared, file.path(res_dir, paste0(PRJNA,"_gsea_shared_terms.csv")), row.names = FALSE)
print(head(shared[, c("Description", "setSize", "NES_field", "NES_gh")], 20))

# ------------------------------------------------------------------
# Over-representation on genes significant in BOTH environments
# ------------------------------------------------------------------
t2g_u <- t2g_exp[t2g_exp$gene %in% universe, ]
run_ora <- function(genes) {
  as.data.frame(enricher(genes, universe = universe, TERM2GENE = t2g_u, TERM2NAME = t2n,
                         minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05))
}
ora_cols <- c("Description", "Count", "GeneRatio", "p.adjust")

ora_up   <- run_ora(shared_up)
ora_down <- run_ora(shared_down)
write.csv(ora_up,   file.path(res_dir, "ora_shared_up.csv"),   row.names = FALSE)
write.csv(ora_down, file.path(res_dir, "ora_shared_down.csv"), row.names = FALSE)

cat("\nShared UP:\n");   print(head(ora_up[, ora_cols], 15))
cat("\nShared DOWN:\n"); print(head(ora_down[, ora_cols], 15))

cat("\nDone. Outputs in", res_dir, "\n")