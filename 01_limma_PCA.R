# Limma-voom differential expression + PCA
# Contrasts: FLvsHL (FL_T30 vs HL_T30), FLvsDark (FL_T30 vs Dark_T30)

library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)

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

# 11. PCA on voom log-CPM values
pca_res <- prcomp(t(v$E), scale. = FALSE)
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df  <- as.data.frame(pca_res$x[, 1:18])
pca_df$condition <- sample_info$group
pca_df$sample    <- sample_info$sample

cond_colors <- c("T0"       = "#606060",
                  "Dark_T30" = "black",
                  "FL_T30"   = "#4C8FCD",
                  "LL_T30"   = "#4DAF4A",
                  "ML_T30"   = "#E8821A",
                  "HL_T30"   = "#CC3333")

cond_labels <- c("T0"       = "T0",
                  "Dark_T30" = "Dark",
                  "FL_T30"   = "FL",
                  "LL_T30"   = "LL",
                  "ML_T30"   = "ML",
                  "HL_T30"   = "HL")

pca_theme <- theme_classic(base_size = 13) +
  theme(legend.position = "right",
        legend.title    = element_text(face = "bold"),
        legend.key.size = unit(0.9, "lines"),
        plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey40"),
        axis.title      = element_text(size = 12),
        panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5))

# e.g. "Dark_T30_1" -> "Rep 1"
make_rep_label <- function(s) {
  num <- sub(".*_([0-9]+)$", "\\1", s)
  paste0("Rep ", num)
}

make_pca_plot <- function(df, xvar, yvar, xlabel, ylabel, title, subtitle = NULL) {
  df$rep_label <- make_rep_label(df$sample)

  ggplot(df, aes_string(x = xvar, y = yvar, color = "condition")) +
    geom_point(size = 5, alpha = 0.95) +
    geom_text_repel(aes(label = rep_label),
                     size = 3.5, fontface = "plain",
                     max.overlaps = 30,
                     show.legend = FALSE,
                     box.padding = 0.35,
                     point.padding = 0.3) +
    scale_color_manual(values = cond_colors, labels = cond_labels, name = "Condition") +
    labs(x = xlabel, y = ylabel, title = title, subtitle = subtitle) +
    pca_theme
}

# 12. PCA scatter plots (PC1-PC2, PC2-PC3, PC3-PC4)
p_pc12 <- make_pca_plot(pca_df, "PC1", "PC2",
                         paste0("PC1 (", pct_var[1], "%)"),
                         paste0("PC2 (", pct_var[2], "%)"),
                         "PCA - PC1 vs PC2")

# Fixed 1:1 aspect ratio so PC1/PC2 are drawn on the same screen
# scale, reflecting that PC2 explains less variance than PC1
p_pc12 <- p_pc12 + coord_fixed(ratio = 1)
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

# 13. Scree plot (top 10 PCs)
n_pcs_available <- min(10, length(pct_var))
scree_df <- data.frame(
  PC     = paste0("PC", 1:n_pcs_available),
  VarExp = pct_var[1:n_pcs_available],
  CumVar = cumsum(pct_var[1:n_pcs_available])
)
scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)

p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = VarExp), fill = "#4393C3", alpha = 0.85) +
  geom_line(aes(y = CumVar, group = 1), color = "#D6604D", linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = CumVar), color = "#D6604D", size = 2.5) +
  labs(x = "Principal Component", y = "Variance Explained (%)",
       title = "Scree Plot", caption = "Dashed line = cumulative variance") +
  pca_theme +
  theme_bw(base_size = 13)

ggsave("Scree_plot.pdf", p_scree, width = 8, height = 6, dpi = 300)

# 14. Combined PCA figure panel
png("PCA_Plots.png", width = 14, height = 12, units = "in", res = 200, bg = "white")
(p_pc12 | p_pc23) / (p_pc34 | p_scree) +
  plot_annotation(title = "Principal Component Analysis - voom log-CPM values",
                   theme = theme(plot.title = element_text(face = "bold", size = 16)))
dev.off()

# 15. PC1 and PC2 loadings (top contributing genes)
n_top_loadings <- 20

loadings <- as.data.frame(pca_res$rotation[, 1:2])
loadings$Gene_ID <- rownames(loadings)

gene_name_map <- raw %>% select(Gene_ID, Gene_Name)
loadings <- loadings %>%
  left_join(gene_name_map, by = "Gene_ID") %>%
  mutate(label = ifelse(is.na(Gene_Name) | Gene_Name == "", Gene_ID, Gene_Name))

make_loadings_plot <- function(df, pc_col, pc_num, pct) {
  top_df <- df %>%
    arrange(desc(abs(.data[[pc_col]]))) %>%
    slice_head(n = n_top_loadings) %>%
    arrange(.data[[pc_col]]) %>%
    mutate(label = factor(label, levels = label))

  ggplot(top_df, aes(x = .data[[pc_col]], y = label, fill = .data[[pc_col]] > 0)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    scale_fill_manual(values = c(`TRUE` = "#CC3333", `FALSE` = "#4C8FCD")) +
    labs(x = paste0("Loading on PC", pc_num, " (", pct, "%)"),
         y = NULL,
         title = paste0("Top ", n_top_loadings, " genes - PC", pc_num, " loadings")) +
    theme_classic(base_size = 12) +
    theme(plot.title  = element_text(face = "bold", size = 13, hjust = 0.5),
          axis.text.y = element_text(size = 8))
}

p_load_pc1 <- make_loadings_plot(loadings, "PC1", 1, pct_var[1])
p_load_pc2 <- make_loadings_plot(loadings, "PC2", 2, pct_var[2])

ggsave("PC1_loadings.pdf", p_load_pc1, width = 7, height = 8, dpi = 300)
ggsave("PC2_loadings.pdf", p_load_pc2, width = 7, height = 8, dpi = 300)

png("PCA_Loadings.png", width = 14, height = 9, units = "in", res = 200, bg = "white")
p_load_pc1 | p_load_pc2
dev.off()

cat("\n Complete \n")