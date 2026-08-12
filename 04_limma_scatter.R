# Limma-voom differential expression analysis
# Dataset: NRaw_Counts_GO_Terms.xlsx (sheet "Dataset")
# Contrasts: FLvsHL (FL_T30 vs HL_T30), FLvsDark (FL_T30 vs Dark_T30), FLvsML (FL_T30 vs ML_T30)
# Plus: FL/HL vs FL/Dark fold-change scatter comparison

library(readxl)
library(edgeR)
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

# 1. Import data
input_file <- "NRaw_Counts_GO_Terms.xlsx"
raw <- as.data.frame(read_excel(input_file, sheet = "Dataset"))

# 2. Build the count matrix
count_cols <- c(
  "T0_1",       "T0_2",       "T0_3",
  "Dark_T30_1", "Dark_T30_2", "Dark_T30_3",
  "LL_T30_1",   "LL_T30_2",   "LL_T30_3",
  "ML_T30_1",   "ML_T30_2",   "ML_T30_3",
  "HL_T30_1",   "HL_T30_2",   "HL_T30_3",
  "FL_T30_1",   "FL_T30_2",   "FL_T30_3"
)

counts <- as.matrix(raw[, count_cols])
rownames(counts) <- raw$Gene_ID

# 3. Sample metadata
sample_info <- data.frame(sample = colnames(counts), stringsAsFactors = FALSE)
sample_info$group <- sub("_[0-9]+$", "", sample_info$sample)
sample_info$group <- factor(
  sample_info$group,
  levels = c("T0", "Dark_T30", "LL_T30", "ML_T30", "HL_T30", "FL_T30")
)

# 4. DGEList + low-count filtering
dge <- DGEList(counts = counts, group = sample_info$group)

keep <- filterByExpr(dge, group = sample_info$group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("Genes retained after filterByExpr: %d / %d\n", sum(keep), length(keep)))

# 5. TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# 6. Design matrix
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

# 7. Voom transformation
v <- voom(dge, design, plot = FALSE)

# 8. Linear model fit
fit <- lmFit(v, design)

# 9. Contrasts
contrast_matrix <- makeContrasts(
  FLvsHL   = FL_T30 - HL_T30,
  FLvsDark = FL_T30 - Dark_T30,
  FLvsML   = FL_T30 - ML_T30,
  levels   = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 10. Extract results
extract_results <- function(fit_obj, coef_name, annotation) {
  tt <- topTable(fit_obj, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "P")
  tt$Gene_ID <- rownames(tt)
  res <- merge(annotation[, c("Gene_ID", "Gene_Name")], tt, by = "Gene_ID", all.y = TRUE)
  res[order(res$P.Value), ]
}

res_FLvsHL   <- extract_results(fit2, "FLvsHL",   raw)
res_FLvsDark <- extract_results(fit2, "FLvsDark", raw)
res_FLvsML   <- extract_results(fit2, "FLvsML",   raw)

# 11. Export results
write.csv(res_FLvsHL,   "limma_FLvsHL_results.csv",   row.names = FALSE)
write.csv(res_FLvsDark, "limma_FLvsDark_results.csv", row.names = FALSE)
write.csv(res_FLvsML,   "limma_FLvsML_results.csv",   row.names = FALSE)

cat(" - limma_FLvsHL_results.csv   (", nrow(res_FLvsHL),   "genes )\n")
cat(" - limma_FLvsDark_results.csv (", nrow(res_FLvsDark), "genes )\n")
cat(" - limma_FLvsML_results.csv   (", nrow(res_FLvsML),   "genes )\n")

# 12. Summary of DEGs (|logFC| >= 1, adj.P.Val < 0.05)
results_list <- list(FLvsHL = res_FLvsHL, FLvsML = res_FLvsML)

cat("\n── Differentially Expressed Genes (|log2FC| >= 1, padj < 0.05) ──\n")
for (nm in names(results_list)) {
  res <- results_list[[nm]]
  sig <- res %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) >= 1)
  up  <- sum(sig$logFC > 0)
  dn  <- sum(sig$logFC < 0)
  cat(sprintf("  %-12s  UP: %4d  |  DOWN: %4d  |  Total: %4d\n",
              nm, up, dn, nrow(sig)))
}

# 13. FL/HL vs FL/Dark fold-change scatter plot
lfc_cut  <- 1
pval_cut <- 0.05

cat_colours <- c(
  "NS"                 = "grey75",
  "Same direction"     = "#D73027",
  "Opposite direction" = "#4575B4"
)

diff_df <- res_FLvsHL %>%
  select(Gene_ID, Gene_Name,
         lfc_HL    = logFC,
         padj_HL   = adj.P.Val) %>%
  inner_join(
    res_FLvsDark %>%
      select(Gene_ID,
             lfc_Dark  = logFC,
             padj_Dark = adj.P.Val),
    by = "Gene_ID"
  ) %>%
  filter(!is.na(lfc_HL), !is.na(lfc_Dark),
         !is.na(padj_HL), !is.na(padj_Dark)) %>%
  mutate(
    sig_HL   = padj_HL < pval_cut & abs(lfc_HL) >= lfc_cut,
    sig_Dark = padj_Dark < pval_cut & abs(lfc_Dark) >= lfc_cut,
    sig_both = sig_HL & sig_Dark,
    same_dir = sign(lfc_HL) == sign(lfc_Dark),
    category = case_when(
      sig_both & same_dir  ~ "Same direction",
      sig_both & !same_dir ~ "Opposite direction",
      TRUE                 ~ "NS"
    )
  )

diff_counts       <- table(diff_df$category)
diff_legend_labels <- setNames(
  paste0(names(diff_counts), "  (n=", diff_counts, ")"),
  names(diff_counts)
)

top_diff <- diff_df %>%
  filter(category != "NS") %>%
  mutate(dist_origin = sqrt(lfc_Dark^2 + lfc_HL^2)) %>%
  group_by(category) %>%
  slice_max(order_by = dist_origin, n = 20) %>%
  ungroup()

write.csv(top_diff,
          file = "Top_DEG_FLHLvsFLDark_limma.csv",
          row.names = FALSE)

diff_plot <- ggplot(diff_df,
                    aes(x = lfc_Dark, y = lfc_HL, colour = category)) +
  geom_point(data = ~ filter(., category == "NS"),
             size = 0.7, alpha = 0.35, stroke = 0) +
  geom_point(data = ~ filter(., category != "NS"),
             size = 1.6, alpha = 0.75, stroke = 0) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
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
    title    = "Fold-Change comparison: FL/HL vs FL/Dark",
    subtitle = sprintf(
      "Significant if padj < %.2f  •  %d genes",
      pval_cut, nrow(diff_df)),
    x        = expression(log[2]~"FC (FL/Dark)"),
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
  filename = "FC_FLonHL_vs_FLonDark_TB.png",
  plot     = diff_plot,
  width    = 7, height = 6.5, dpi = 300
)
cat("Saved: FC_FLonHL_vs_FLonDark_TB.png\n")

cat("\n── FL-HL vs FL-Dark scatter — gene categories ──\n")
print(diff_counts)

cat("\nAnalyse complete.\n")
