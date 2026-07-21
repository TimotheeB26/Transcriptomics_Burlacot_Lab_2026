# =============================================================================
# DESeq2 + Venn / Euler Diagram Pipeline
# =============================================================================

library(DESeq2)
library(dplyr)
library(tibble)
library(eulerr)
library(ggVennDiagram)
library(RColorBrewer)

# ── Thresholds ────────────────────────────────────────────────
LFC_THRESH  <- 1
PADJ_THRESH <- 0.05

# =============================================================================
# 1. Load raw data (adjust path as needed)
# =============================================================================
# raw <- read.csv("your_counts_file.csv")  # must contain Gene_ID, Gene_Name, and *_Raw.Read.Count columns

# ── 2. Build count matrix ─────────────────────────────────────
count_cols <- grep("_Raw.Read.Count$", names(raw), value = TRUE)

counts_mat <- raw[, count_cols] %>%
  mutate(across(everything(), ~ as.integer(round(.)))) %>%
  as.data.frame()
rownames(counts_mat) <- raw$Gene_ID

# ── 3. Build sample metadata ──────────────────────────────────
sample_names <- count_cols
condition <- sub("WT-(FL|HL|ML|T0)-T?30?-[123]_Raw.Read.Count", "\\1", sample_names)
condition[grepl("T0", sample_names)] <- "T0"

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(condition, levels = c("T0", "FL", "HL", "ML")),
  row.names = sample_names
)
cat("Sample table:\n"); print(col_data)

# ── 4. Create DESeq2 object & run analysis ────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ condition
)

# Pre-filter: keep genes with ≥ 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message(sprintf("\nGènes après filtrage : %d", nrow(dds)))

dds <- DESeq(dds)

# ── 5. Extract results for all three contrasts ────────────────
contrasts <- list(
  FL_vs_T0 = c("condition", "FL", "T0"),
  HL_vs_T0 = c("condition", "HL", "T0"),
  ML_vs_T0 = c("condition", "ML", "T0")
)

results_list <- lapply(names(contrasts), function(nm) {
  res <- results(dds,
                 contrast      = contrasts[[nm]],
                 alpha         = 0.05,
                 pAdjustMethod = "BH") %>%
    as.data.frame() %>%
    rownames_to_column("Gene_ID") %>%
    left_join(raw[, c("Gene_ID", "Gene_Name")], by = "Gene_ID") %>%
    mutate(comparison = nm)
  res
})
names(results_list) <- names(contrasts)

# ── 6. Build DEG lists (|LFC| ≥ 1, FDR < 0.05) ───────────────
deg_list <- lapply(results_list, function(res) {
  res %>%
    filter(!is.na(padj),
           padj < PADJ_THRESH,
           abs(log2FoldChange) >= LFC_THRESH) %>%
    pull(Gene_ID)
})
names(deg_list) <- c("FL vs T0", "HL vs T0", "ML vs T0")

cat("\nDEG counts per contrast:\n")
print(sapply(deg_list, length))

# =============================================================================
# 6. Proportional Venn / Euler Diagram (eulerr)
# =============================================================================
cat("\nGenerating proportional Euler diagram...\n")

euler_fit <- euler(deg_list, shape = "ellipse")

venn_cols <- brewer.pal(3, "Set2")

png("Venn_DEGs_proportional.png", width = 2400, height = 2000, res = 300)
plot(euler_fit,
     fills      = list(fill = venn_cols, alpha = 0.55),
     edges      = list(col = venn_cols, lwd = 2),
     labels     = list(font = 2, cex = 1.2,
                       labels = c("FL vs T0", "HL vs T0", "ML vs T0")),
     quantities = list(font = 1, cex = 1.0, type = "counts"),
     main       = list(label = "Venn diagram - DEG per condition vs T0  (|LFC| ≥ 1, padj < 0.05)",
                       font = 2, cex = 1.3))
dev.off()
cat("Proportional Venn diagram saved: Venn_DEGs_proportional.png\n")

# =============================================================================
# 7. Classic Venn with ggVennDiagram
# =============================================================================
cat("Generating standard Venn diagram (ggVennDiagram)...\n")

p_venn <- ggVennDiagram(
  deg_list,
  label_alpha = 0,
  label       = "count",
  set_color   = venn_cols
) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title    = "DEGs per Condition vs T0",
       subtitle = sprintf("|LFC| ≥ %.0f, FDR < %.2f", LFC_THRESH, PADJ_THRESH),
       fill     = "Gene count") +
  theme_classic() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("Venn_DEGs_ggVennDiagram.png", p_venn,
       width = 8, height = 7, dpi = 300)
cat("Standard Venn diagram saved: Venn_DEGs_ggVennDiagram.png\n")

# =============================================================================
# 10. Intersection summary table & per-region gene lists
# =============================================================================
cat("\n--- Intersection summary ---\n")

FL <- deg_list[["FL vs T0"]]
HL <- deg_list[["HL vs T0"]]
ML <- deg_list[["ML vs T0"]]

summary_df <- data.frame(
  Region     = c("FL only", "HL only", "ML only",
                 "FL ∩ HL", "FL ∩ ML", "HL ∩ ML",
                 "FL ∩ HL ∩ ML"),
  Gene_count = c(
    length(setdiff(FL, union(HL, ML))),
    length(setdiff(HL, union(FL, ML))),
    length(setdiff(ML, union(FL, HL))),
    length(setdiff(intersect(FL, HL), ML)),
    length(setdiff(intersect(FL, ML), HL)),
    length(setdiff(intersect(HL, ML), FL)),
    length(Reduce(intersect, list(FL, HL, ML)))
  )
)

print(summary_df)
write.csv(summary_df, "DEG_intersection_summary.csv", row.names = FALSE)

# Save gene lists per intersection region
intersect_genes <- list(
  FL_only  = setdiff(FL, union(HL, ML)),
  HL_only  = setdiff(HL, union(FL, ML)),
  ML_only  = setdiff(ML, union(FL, HL)),
  FL_HL    = setdiff(intersect(FL, HL), ML),
  FL_ML    = setdiff(intersect(FL, ML), HL),
  HL_ML    = setdiff(intersect(HL, ML), FL),
  FL_HL_ML = Reduce(intersect, list(FL, HL, ML))
)

dir.create("DESeq2_results", showWarnings = FALSE)

gene_annotation <- raw[, c("Gene_ID", "Gene_Name")]

for (nm in names(intersect_genes)) {
  ids <- intersect_genes[[nm]]
  if (length(ids) > 0) {
    ann <- gene_annotation %>% filter(Gene_ID %in% ids)
    write.csv(ann,
              file      = sprintf("DESeq2_results/Genes_%s.csv", nm),
              row.names = FALSE)
  }
}

cat("\nAll outputs saved.\n")
cat("Key files:\n")
cat("  Venn_DEGs_proportional.png   – proportional Euler diagram (eulerr)\n")
cat("  Venn_DEGs_ggVennDiagram.png  – standard Venn (ggVennDiagram)\n")
cat("  DEG_intersection_summary.csv – count per Venn region\n")
cat("  DESeq2_results/              – gene lists per intersection\n")