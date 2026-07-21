#
#  DESeq2 Analysis + Volcano Plots
#  Dataset: WT light conditions (Dark, FL, HL, LL, ML at T30; HL also has T0) - HTSeq counts
#  Contrasts of interest: FL vs HL (both at T30) and FL vs Dark (both at T30)
#

# Install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("DESeq2", "ggplot2", "ggrepel", "readxl", "dplyr", "patchwork", "ggnewscale")
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
library(ggnewscale)

# Run if 'select' error
# install.packages("conflicted")
# library(conflicted)
# conflict_prefer("select", "dplyr")
# conflict_prefer("filter", "dplyr")


# 1. Load dataset

data_path <- "Raw_Normalized_Counts_Genes_DESeq2_HTSEQ_New.xlsx"
raw <- read_excel(data_path, sheet = "All.Genes")

# Round to integers for DESeq2

count_cols <- grep("_Raw.Read.Count$", names(raw), value = TRUE)
counts_mat <- raw[, count_cols] %>%
  mutate(across(everything(), ~ as.integer(round(.)))) %>%
  as.data.frame()

rownames(counts_mat) <- raw$Gene_ID


# 2. Build sample metadata
#    Sample names look like: Dark-T30-1_Raw.Read.Count, FL-T30-2_Raw.Read.Count,
#    HL-T0-1_Raw.Read.Count, HL-T30-3_Raw.Read.Count, LL-T30-1_..., ML-T30-1_...
#    NOTE: HL appears at both T0 and T30, so condition must encode light AND
#    timepoint (otherwise HL-T0 and HL-T30 samples would be merged together).

sample_names <- count_cols

light <- sub("^(Dark|FL|HL|LL|ML)-T(0|30)-[123]_Raw.Read.Count$", "\\1", sample_names)
time  <- sub("^(Dark|FL|HL|LL|ML)-T(0|30)-[123]_Raw.Read.Count$", "T\\2", sample_names)

condition <- paste(light, time, sep = "_")

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(condition, levels = c("HL_T0", "Dark_T30", "FL_T30", "HL_T30", "LL_T30", "ML_T30")),
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


# 4. Extract results for the two contrasts of interest (both at T30)

contrasts <- list(
  FL_vs_HL   = c("condition", "FL_T30", "HL_T30"),
  FL_vs_Dark = c("condition", "FL_T30", "Dark_T30")
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


# 4b. Genes of interest to highlight on every volcano plot
#     (named vector: Gene_ID = display label)

highlight_genes <- c(
  "Cre12.g531900" = "FLVA",
  "Cre16.g661200" = "THB8",
  "Cre13.g588150" = "VTC2",
  "Cre15.g635800" = "SMC1"
)


# 5. Volcano plot function

make_volcano <- function(res_df, title, lfc_cut  = 1, pval_cut = 0.05, top_n = 15,
                          highlight = highlight_genes) {
  df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10p = -log10(pmax(padj, 1e-300)),
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

  # Genes of interest present in this comparison's results
  # NOTE: Cre-locus IDs (e.g. "Cre12.g531900") live in the Gene_Name column,
  # not Gene_ID (which holds the CHLRE_v5 identifier)
  highlight_df <- df %>%
    filter(Gene_Name %in% names(highlight)) %>%
    mutate(label = highlight[Gene_Name])

  missing_genes <- setdiff(names(highlight), df$Gene_Name)
  if (length(missing_genes) > 0) {
    warning(sprintf(
      "[%s] Highlight gene(s) not found in tested results (filtered out or absent): %s",
      title, paste(missing_genes, collapse = ", ")
    ))
  }

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10p, colour = significance)) +
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
               linetype = "dashed", colour = "grey40", linewidth = 0.5)

  # Overlay genes of interest: distinct shape/colour + repelled labels,
  # drawn on top of and independently from the significance colour scale
  if (nrow(highlight_df) > 0) {
    p <- p +
      ggnewscale::new_scale_colour() +
      geom_point(
        data = highlight_df,
        aes(x = log2FoldChange, y = neg_log10p),
        shape = 21, fill = "yellow", colour = "black",
        size = 3.2, stroke = 0.8, inherit.aes = FALSE
      ) +
      geom_text_repel(
        data = highlight_df,
        aes(x = log2FoldChange, y = neg_log10p, label = label),
        inherit.aes = FALSE,
        fontface = "bold", size = 3.6, colour = "black",
        box.padding = 0.5, point.padding = 0.3,
        segment.colour = "black", segment.size = 0.4,
        max.overlaps = Inf, min.segment.length = 0
      )
  }

  p +
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
  FL_vs_HL   = "FL vs HL",
  FL_vs_Dark = "FL vs Dark"
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

combined <- wrap_plots(volcano_plots, ncol = 2) +
  plot_annotation(
    title    = "DESeq2 Volcano Plots - FL vs HL and FL vs Dark",
    subtitle = "Thresholds: |log2FC| >= 1 AND padj < 0.05",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40")
    )
  )

ggsave(
  filename = "Volcano_AllComparisons_Panel.png",
  plot     = combined,
  width    = 14, height = 7, dpi = 300
)

cat("Analysis complete")
