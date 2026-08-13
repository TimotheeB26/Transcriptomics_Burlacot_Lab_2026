# Limma-voom differential expression analysis + volcano plots
# Contrasts: DarkvsT0, LLvsT0, MLvsT0, HLvsDark, FLvsT0, FLvsHL, FLvsDark

library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(patchwork)
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
  DarkvsT0 = Dark_T30 - T0,
  LLvsT0   = LL_T30   - T0,
  MLvsT0   = ML_T30   - T0,
  HLvsDark = HL_T30   - Dark_T30,
  FLvsT0   = FL_T30   - T0,
  FLvsHL   = FL_T30   - HL_T30,
  FLvsDark = FL_T30   - Dark_T30,
  levels   = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 10. Extract results
extract_results <- function(fit_obj, coef_name, annotation) {
  tt <- topTable(fit_obj, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "P")
  tt$Gene_ID <- rownames(tt)
  res <- merge(annotation[, c("Gene_ID", "Gene_Name", "Gene_Symbol")], tt, by = "Gene_ID", all.y = TRUE)
  res[order(res$P.Value), ]
}

res_DarkvsT0 <- extract_results(fit2, "DarkvsT0", raw)
res_LLvsT0   <- extract_results(fit2, "LLvsT0",   raw)
res_MLvsT0   <- extract_results(fit2, "MLvsT0",   raw)
res_HLvsDark <- extract_results(fit2, "HLvsDark", raw)
res_FLvsT0   <- extract_results(fit2, "FLvsT0",   raw)
res_FLvsHL   <- extract_results(fit2, "FLvsHL",   raw)
res_FLvsDark <- extract_results(fit2, "FLvsDark", raw)

# 11. Export results
# All outputs (CSVs + volcano plots) are grouped into a single output folder.
output_dir <- "limma_volcano_outputs"
dir.create(output_dir, showWarnings = FALSE)

write.csv(res_DarkvsT0, file.path(output_dir, "limma_DarkvsT0_results.csv"), row.names = FALSE)
write.csv(res_LLvsT0,   file.path(output_dir, "limma_LLvsT0_results.csv"),   row.names = FALSE)
write.csv(res_MLvsT0,   file.path(output_dir, "limma_MLvsT0_results.csv"),   row.names = FALSE)
write.csv(res_HLvsDark, file.path(output_dir, "limma_HLvsDark_results.csv"), row.names = FALSE)
write.csv(res_FLvsT0,   file.path(output_dir, "limma_FLvsT0_results.csv"),   row.names = FALSE)
write.csv(res_FLvsHL,   file.path(output_dir, "limma_FLvsHL_results.csv"),   row.names = FALSE)
write.csv(res_FLvsDark, file.path(output_dir, "limma_FLvsDark_results.csv"), row.names = FALSE)

# 12. Volcano plot builder
# Both the y-axis and the significance call use the BH-adjusted p-value
# (adj.P.Val), not the raw P.Value, so what's plotted matches what's
# used to threshold significance.
make_volcano <- function(res_df, title, contrast_name = title, lfc_cut = 1, pval_cut = 0.05) {
  df <- res_df %>%
    filter(!is.na(logFC)) %>%
    mutate(
      neg_log10p = -log10(pmax(adj.P.Val, 1e-300)),
      significance = case_when(
        adj.P.Val < pval_cut & logFC >=  lfc_cut ~ "Up",
        adj.P.Val < pval_cut & logFC <= -lfc_cut ~ "Down",
        TRUE                                      ~ "NS"
      )
    )

  n_up   <- sum(df$significance == "Up")
  n_down <- sum(df$significance == "Down")

  # Label priority: Gene_Symbol first, falling back to Gene_Name, then Gene_ID
  has_symbol <- "Gene_Symbol" %in% names(df)
  has_name   <- "Gene_Name"   %in% names(df)

  sym  <- if (has_symbol) df$Gene_Symbol else NA_character_
  name <- if (has_name)   df$Gene_Name   else NA_character_

  df$label <- dplyr::case_when(
    !is.na(sym)  & sym  != "" ~ sym,
    !is.na(name) & name != "" ~ name,
    TRUE                      ~ df$Gene_ID
  )

  # Top 10 up- and down-regulated genes, ranked by adj.P.Val
  top_up <- df %>%
    filter(significance == "Up") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 10)

  top_down <- df %>%
    filter(significance == "Down") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 10)

  top_labels <- bind_rows(top_up, top_down)

  # Export the genes labeled on the plot (Gene_Symbol prioritized, Gene_Name fallback)
  if (nrow(top_labels) > 0) {
    label_export <- top_labels[, c("Gene_ID", "label", "logFC", "adj.P.Val", "significance")]
    names(label_export)[names(label_export) == "label"] <- "Labeled_Gene"
    write.csv(
      label_export,
      file.path(output_dir, paste0("Volcano_", contrast_name, "_labeled_genes.csv")),
      row.names = FALSE
    )
  }

  ggplot(df, aes(x = logFC, y = neg_log10p, colour = significance)) +
    geom_point(aes(size = significance), alpha = 0.6, stroke = 0) +
    scale_size_manual(values = c(Up = 1.8, Down = 1.8, NS = 0.8)) +
    scale_colour_manual(
      values = c(Up = "#E63946", Down = "#457B9D", NS = "grey70"),
      labels = c(
        Up   = paste0("Up (", n_up, ")"),
        Down = paste0("Down (", n_down, ")"),
        NS   = "Not significant"
      )
    ) +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text_repel(
      data = top_labels,
      aes(label = label),
      size = 2.8,
      colour = "black",
      fontface = "italic",
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.3,
      segment.colour = "grey30",
      box.padding = 0.35,
      point.padding = 0.2,
      seed = 42,
      show.legend = FALSE
    ) +
    labs(
      title    = title,
      subtitle = sprintf("|log2FC| >= %g, padj < %g  -  %d genes tested", lfc_cut, pval_cut, nrow(df)),
      x        = expression(log[2]~("Fold Change")),
      y        = expression(-log[10]~("padj")),
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
    guides(size = "none", colour = guide_legend(override.aes = list(size = 3)))
}

# 13. Build and save individual volcano plots
results_list <- list(
  DarkvsT0 = res_DarkvsT0,
  LLvsT0   = res_LLvsT0,
  MLvsT0   = res_MLvsT0,
  HLvsDark = res_HLvsDark,
  FLvsT0   = res_FLvsT0,
  FLvsHL   = res_FLvsHL,
  FLvsDark = res_FLvsDark
)

plot_titles <- c(
  DarkvsT0 = "Dark vs T0",
  LLvsT0   = "LL vs T0",
  MLvsT0   = "ML vs T0",
  HLvsDark = "HL vs Dark",
  FLvsT0   = "FL vs T0",
  FLvsHL   = "FL vs HL",
  FLvsDark = "FL vs Dark"
)

volcano_plots <- lapply(names(results_list), function(nm) {
  make_volcano(results_list[[nm]], title = plot_titles[nm], contrast_name = nm)
})
names(volcano_plots) <- names(results_list)

for (nm in names(volcano_plots)) {
  ggsave(
    filename = file.path(output_dir, paste0("Volcano_", nm, ".png")),
    plot     = volcano_plots[[nm]],
    width    = 7, height = 6, dpi = 300
  )
}

# 14. Combined panel of all contrasts
combined <- wrap_plots(volcano_plots, ncol = 4) +
  plot_annotation(
    title    = "limma-voom Volcano Plots - All Contrasts",
    subtitle = "Thresholds: |log2FC| >= 1 AND padj < 0.05",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40")
    )
  )

ggsave(
  filename = file.path(output_dir, "Volcano_AllComparisons_Panel.png"),
  plot     = combined,
  width    = 26, height = 12, dpi = 300
)

cat("\n Complete \n")