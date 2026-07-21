## ============================================================
##  DESeq2 Transcriptomic Analysis Pipeline
##  Conditions: WT-FL-T30, WT-HL-T30, WT-ML-T30, WT-T0
## ============================================================

## ---- 0. Install / load packages ----------------------------
pkgs <- c("DESeq2", "ggplot2", "ggrepel",
          "RColorBrewer", "readxl", "dplyr", "patchwork")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("DESeq2")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(p)
    } else {
      install.packages(p)
    }
  }
  library(p, character.only = TRUE)
}


## ---- 1. Load data ------------------------------------------
raw <- read_excel("Raw_Normalized_Counts_Genes_DESeq2_ToUse.xlsx",
                  sheet = "All.Genes")

# Gene IDs as row names; keep only integer count columns
count_mat <- raw %>%
  column_to_rownames("Gene_ID") %>%
  select(-Gene_Name) %>%
  as.matrix()

storage.mode(count_mat) <- "integer"

# Remove genes with zero counts in all samples
count_mat <- count_mat[rowSums(count_mat) > 0, ]
cat("Genes after removing all-zero rows:", nrow(count_mat), "\n")


## ---- 2. Build colData (sample metadata) --------------------
sample_names <- colnames(count_mat)

# Parse condition from column names  (strip _Raw.Read.Count)
condition <- sub("_Raw\\.Read\\.Count", "", sample_names)
condition <- sub("-[0-9]+$", "", condition)          # drop replicate number

colData <- data.frame(
  row.names   = sample_names,
  condition   = factor(condition,
                        levels = c("WT-T0", "WT-FL-T30",
                                   "WT-HL-T30", "WT-ML-T30"))
)
cat("\nSample table:\n")
print(colData)


## ---- 3. DESeq2 object & pre-filtering ----------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = colData,
  design    = ~ condition
)

# Keep genes with >= 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("\nGenes after low-count filtering:", nrow(dds), "\n")


## ---- 4. Run DESeq2 -----------------------------------------
dds <- DESeq(dds)


## ---- 5. Variance-stabilising transformation for QC plots ---
vst_data <- vst(dds, blind = TRUE)


## ============================================================
##  FIGURE 1 – PCA (PC1 vs PC2, PC2 vs PC3, PC3 vs PC4, Scree)
## ============================================================

pca_res  <- prcomp(t(assay(vst_data)), scale. = FALSE)
pct_var  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df   <- as.data.frame(pca_res$x[, 1:4])
pca_df$condition <- colData$condition
pca_df$sample    <- rownames(pca_df)

cond_colors <- c("WT-T0"     = "#606060",   # dark grey
                 "WT-FL-T30" = "#4C8FCD",   # steel blue
                 "WT-HL-T30" = "#CC3333",   # red
                 "WT-ML-T30" = "#E8821A")   # orange

# Short condition labels for the legend
cond_labels <- c("WT-T0"     = "T0",
                 "WT-FL-T30" = "FL",
                 "WT-HL-T30" = "HL",
                 "WT-ML-T30" = "ML")

pca_theme <- theme_classic(base_size = 13) +
  theme(legend.position      = "right",
        legend.title         = element_text(face = "bold"),
        legend.key.size      = unit(0.9, "lines"),
        plot.title           = element_text(face = "bold", size = 14,
                                            hjust = 0.5),
        plot.subtitle        = element_text(size = 9, hjust = 0.5,
                                            color = "grey40"),
        axis.title           = element_text(size = 12),
        panel.border         = element_rect(color = "grey70",
                                            fill = NA, linewidth = 0.5))

# Rep labels: strip condition prefix, keep "Rep N"
make_rep_label <- function(s) {
  # e.g. "WT-FL-T30-1_Raw.Read.Count" -> "Rep 1"
  num <- sub(".*-(\\d+)_Raw\\.Read\\.Count", "\\1", s)
  paste0("Rep ", num)
}

make_pca_plot <- function(df, xvar, yvar, xlabel, ylabel, title, subtitle = NULL) {
  df$rep_label <- make_rep_label(df$sample)

  ggplot(df, aes_string(x = xvar, y = yvar, color = "condition")) +
    geom_point(size = 5, alpha = 0.95) +
    geom_text_repel(aes(label = rep_label),
                    size = 3.5, fontface = "plain",
                    max.overlaps = 30,
                    show.legend  = FALSE,
                    box.padding  = 0.35,
                    point.padding = 0.3) +
    scale_color_manual(values = cond_colors,
                       labels = cond_labels,
                       name   = "Condition") +
    labs(x = xlabel, y = ylabel,
         title    = title,
         subtitle = subtitle) +
    pca_theme
}


p_pc12 <- make_pca_plot(pca_df, "PC1", "PC2",
  paste0("PC1 (", pct_var[1], "%)"),
  paste0("PC2 (", pct_var[2], "%)"),
  "PCA - PC1 vs PC2")

ggsave("PC1vsPC2.pdf", p_pc12, width = 8, height = 6, dpi = 300)

p_pc23 <- make_pca_plot(pca_df, "PC2", "PC3",
  paste0("PC2 (", pct_var[2], "%)"),
  paste0("PC3 (", pct_var[3], "%)"),
  "PCA - PC2 vs PC3")

ggsave("PC2vsPC3.pdf", p_pc23, width = 8, height = 6, dpi = 300)

p_pc34 <- make_pca_plot(pca_df, "PC3", "PC4",
                        paste0("PC3 (", pct_var[3], "%)"),
                        paste0("PC4 (", pct_var[4], "%)"),
                        "PCA - PC3 vs PC4")


ggsave("PC3vsPC4.pdf", p_pc34, width = 8, height = 6, dpi = 300)


# Scree plot (top 10 PCs)
scree_df <- data.frame(
  PC      = paste0("PC", 1:10),
  VarExp  = pct_var[1:10],
  CumVar  = cumsum(pct_var[1:10])
)
scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)

p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = VarExp), fill = "#4393C3", alpha = 0.85) +
  geom_line(aes(y = CumVar, group = 1), color = "#D6604D",
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = CumVar), color = "#D6604D", size = 2.5) +
  labs(x = "Principal Component",
       y = "Variance Explained (%)",
       title  = "Scree Plot",
       caption = "Dashed line = cumulative variance") +
  pca_theme +
  theme_bw(base_size = 13)

ggsave("Scree_plot.pdf", p_scree, width = 8, height = 6, dpi = 300)


png("PCA_Plots.png", width = 14, height = 12, units = "in",
    res = 200, bg = "white")
(p_pc12 | p_pc23) / (p_pc34 | p_scree) +
  plot_annotation(title = "Principal Component Analysis – VST-transformed counts",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))
dev.off()
cat("Saved: PCA_Plots.png\n")