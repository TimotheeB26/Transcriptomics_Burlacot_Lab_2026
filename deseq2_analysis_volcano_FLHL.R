#
#  DESeq2 Analysis + Volcano Plots
#  Dataset: WT light conditions (FL, HL, ML at T30) vs T0 (controlled conditions)
#


# Install packages

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

# Run if 'select' error
# install.packages("conflicted")
# library(conflicted)
# conflict_prefer("select", "dplyr")
# conflict_prefer("filter", "dplyr")


# 1. Load dataset

data_path <- "Raw_Normalized_Counts_Genes_DESeq2_ToUse.xlsx"
raw <- read_excel(data_path, sheet = "All.Genes")

# Round to integers for DESeq2

count_cols <- grep("_Raw.Read.Count$", names(raw), value = TRUE)
counts_mat <- raw[, count_cols] %>%
  mutate(across(everything(), ~ as.integer(round(.)))) %>%
  as.data.frame()

rownames(counts_mat) <- raw$Gene_ID


# 2. Build sample metadata

sample_names <- count_cols
condition <- sub("WT-(FL|HL|ML|T0)-T?30?-[123]_Raw.Read.Count", "\\1", sample_names)
condition[grepl("T0", sample_names)] <- "T0"

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(condition, levels = c("T0", "FL", "HL", "ML")),
  row.names = sample_names
)


# 3. Create DESeq2 object and run analysis

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ condition
)

# Manual filter: keep genes with >= 10 counts in at least 3 samples

keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

dds <- DESeq(dds)


# 4. Extract results for all three contrasts

contrasts <- list(
  FL_vs_HL = c("condition", "FL", "HL")
)

results_list <- lapply(names(contrasts), function(nm) {
  res <- results(dds, contrast  = contrasts[[nm]], alpha = 0.05, pAdjustMethod = "BH") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene_ID") %>%
    left_join(raw[, c("Gene_ID", "Gene_Name")], by = "Gene_ID") %>%
    mutate(comparison = nm)
  res
})
names(results_list) <- names(contrasts)

# Save results to CSV files

for (nm in names(results_list)) {
  write.csv(results_list[[nm]],
            file = paste0("DESeq2_results_", nm, ".csv"),
            row.names = FALSE)
}


# 5. Volcano plot function

make_volcano <- function(res_df, title, lfc_cut  = 1, pval_cut = 0.05, top_n = 15) {
  df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10p = -log10(pmax(pvalue, 1e-300)),
      significance = case_when(
        padj < pval_cut & log2FoldChange >=  lfc_cut ~ "Up",
        padj < pval_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE                                          ~ "NS"
      )
    )
  
  # Top genes to label (most significant among DEGs)
  top_genes <- df %>%
    filter(significance != "NS") %>%
    arrange(padj) %>%
    slice_head(n = top_n)
  
  n_up   <- sum(df$significance == "Up")
  n_down <- sum(df$significance == "Down")
  
  ggplot(df, aes(x = log2FoldChange, y = neg_log10p, colour = significance)) +
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
    geom_hline(yintercept = -log10(pval_cut),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    labs(
      title    = title,
      subtitle = sprintf("|log2FC| >= %g, padj < %g  -  %d genes tested",
                         lfc_cut, pval_cut, nrow(df)),
      x        = expression(log[2]~("Fold Change")),
      y        = expression(-log[10]~(padj)),
      colour   = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "top",
      legend.text     = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(size = "none", colour = guide_legend(override.aes = list(size = 3)))
}

# 6. Generate and save individual volcano plots

plot_titles <- c(
  FL_vs_HL = "FL vs HL"
)

volcano_plots <- lapply(names(results_list), function(nm) {
  make_volcano(results_list[[nm]], title = plot_titles[nm])
})
names(volcano_plots) <- names(results_list)

# Save each individual plot

for (nm in names(volcano_plots)) {
  ggsave(
    filename = paste0("Volcano_", nm, ".png"),
    plot     = volcano_plots[[nm]],
    width    = 7, height = 6, dpi = 300
  )
}
# 7. Combine figures

combined <- wrap_plots(volcano_plots, ncol = 3) +
  plot_annotation(
    title    = "DESeq2 Volcano Plots - FL vs HL",
    subtitle = "Thresholds: |log2FC| >= 1 AND padj < 0.05",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40")
    )
  )

ggsave(
  filename = "Volcano_AllComparisons_Panel.png",
  plot     = combined,
  width    = 20, height = 7, dpi = 300
)

cat("Analysis complete")