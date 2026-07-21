# ============================================================
#  DESeq2 Differential Expression Analysis + Volcano Plots
#  Dataset: WT conditions (FL, HL, ML at T30) vs T0 (control)
# ============================================================

# ── 0. Install / load packages ────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("DESeq2", "ggplot2", "ggrepel", "readxl", "dplyr", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("DESeq2")) BiocManager::install(p, ask = FALSE)
    else install.packages(p)
  }
}

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)
library(patchwork)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ── 1. Load data ──────────────────────────────────────────────
data_path <- "Raw_Normalized_Counts_Genes_DESeq2_ToUse.xlsx"   # adjust if needed

raw <- read_excel(data_path, sheet = "All.Genes")

# Round to integers (DESeq2 requires integer counts)
count_cols <- grep("_Raw.Read.Count$", names(raw), value = TRUE)
counts_mat <- raw[, count_cols] %>%
  mutate(across(everything(), ~ as.integer(round(.)))) %>%
  as.data.frame()

rownames(counts_mat) <- raw$Gene_ID

# ── 2. Build sample metadata ──────────────────────────────────
# Conditions inferred from column names
sample_names <- count_cols
condition <- sub("WT-(FL|HL|ML|T0)-T?30?-[123]_Raw.Read.Count", "\\1", sample_names)
condition[grepl("T0", sample_names)] <- "T0"

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(condition, levels = c("T0", "FL", "HL", "ML")),
  row.names = sample_names
)

cat("Sample table:\n")
print(col_data)

# ── 3. Create DESeq2 object & run analysis ────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ condition
)

# Pre-filter: keep genes with ≥ 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat(sprintf("\nGenes after filtering: %d\n", nrow(dds)))

dds <- DESeq(dds)





# ── 4bis. Extract results for all three contrasts ────────────────
contrasts <- list(
  FL_vs_HL = c("condition", "FL", "HL"),
  FL_vs_ML = c("condition", "FL", "ML")
)

results_list <- lapply(names(contrasts), function(nm) {
  res <- results(dds,
                 contrast  = contrasts[[nm]],
                 alpha     = 0.05,
                 pAdjustMethod = "BH") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene_ID") %>%
    left_join(raw[, c("Gene_ID", "Gene_Name")], by = "Gene_ID") %>%
    mutate(comparison = nm)
  res
})
names(results_list) <- names(contrasts)

# Save full results to CSV files
for (nm in names(results_list)) {
  write.csv(results_list[[nm]],
            file = paste0("DESeq2_results_", nm, ".csv"),
            row.names = FALSE)
  cat(sprintf("Saved: DESeq2_results_%s.csv\n", nm))
}

# ── 5. Summary of DEGs ────────────────────────────────────────
cat("\n── Differentially Expressed Genes (|log2FC| ≥ 1, padj < 0.05) ──\n")
for (nm in names(results_list)) {
  res <- results_list[[nm]]
  sig <- res %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  up  <- sum(sig$log2FoldChange > 0)
  dn  <- sum(sig$log2FoldChange < 0)
  cat(sprintf("  %-12s  UP: %4d  |  DOWN: %4d  |  Total: %4d\n",
              nm, up, dn, nrow(sig)))
}




# ── 10. FL/HL vs FL/ML fold-change difference scatter ─────────
# Compute log2FC differences:
#   x = log2FC(FL/HL)
#   y = log2FC(FL/ML)

lfc_cut  <- 1
pval_cut <- 0.05

cat_colours <- c(
  "NS"             = "grey75",
  "Depleted in FL" = "#4575B4",
  "Enriched in FL" = "#D73027"
)


diff_df <- results_list[["FL_vs_HL"]] %>%
  select(Gene_ID, Gene_Name,
         lfc_HL  = log2FoldChange,
         padj_HL = padj) %>%
  inner_join(
    results_list[["FL_vs_ML"]] %>%
      select(Gene_ID, lfc_ML = log2FoldChange, padj_ML = padj),
    by = "Gene_ID"
  ) %>%
  filter(!is.na(lfc_HL), !is.na(lfc_ML),
         !is.na(padj_HL), !is.na(padj_ML)) %>%
  mutate(
    sig_HL = padj_HL < pval_cut,
    sig_ML = padj_ML < pval_cut,
    category = case_when(
      sig_HL & sig_ML & lfc_HL < 0 & lfc_ML < 0 ~ "Depleted in FL",
      sig_HL & sig_ML & lfc_HL > 0 & lfc_ML > 0 ~ "Enriched in FL",
      TRUE                                        ~ "NS"
    )
  )

# Count per category
diff_counts <- table(diff_df$category)
diff_legend_labels <- setNames(
  paste0(names(diff_counts), "  (n=", diff_counts, ")"),
  names(diff_counts)
)

# Top 20 opposite-direction genes furthest from origin, 10 per side
top_diff <- diff_df %>%
  filter(category == "Depleted in FL" | category == "Enriched in FL") %>%
  mutate(dist_origin = sqrt(lfc_ML^2 + lfc_HL^2),
         side = ifelse(lfc_ML > 0, "FL_ML_pos", "FL_ML_neg")) %>%
  group_by(side) %>%
  slice_max(order_by = dist_origin, n = 20) %>%
  ungroup()

write.csv(top_diff,
          file = paste0("Top_DEG_FLHLvsFLML_DESeq.csv"),
          row.names = FALSE)


diff_plot <- ggplot(diff_df,
                    aes(x = lfc_ML, y = lfc_HL, colour = category)) +
  geom_point(data = ~ filter(., category == "NS"),
             size = 0.7, alpha = 0.35, stroke = 0) +
  geom_point(data = ~ filter(., category != "NS"),
             size = 1.6, alpha = 0.75, stroke = 0) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_text_repel(
    data           = top_diff,
    aes(label      = ifelse(!is.na(Gene_Name) & Gene_Name != "",
                            Gene_Name, Gene_ID)),
    size           = 2.6,
    colour         = "black",
    max.overlaps   = 25,
    segment.size   = 0.3,
    segment.colour = "grey50",
    box.padding    = 0.4,
    point.padding  = 0.3
  ) +
  scale_colour_manual(values = cat_colours, labels = diff_legend_labels) +
  labs(
    title    = "Fold-Change comparison: FL/HL vs FL/ML",
    subtitle = sprintf(
      "Significant if padj < %.2f -- %d genes",
      pval_cut, nrow(diff_df)),
    x        = expression(log[2]~"FC (FL/ML)"),
    y        = expression(log[2]~"FC (FL/HL)"),
    colour   = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(
  filename = "FC_FLonHL_vs_FLonML_TB.png",
  plot     = diff_plot,
  width    = 7, height = 6.5, dpi = 300
)
cat("Saved: FC_FLonHL_vs_FLonML_TB.png\n")

cat("\n── FL-HL vs FL-ML scatter — gene categories ──\n")
print(diff_counts)

cat("\n✓ Analysis complete.\n")
