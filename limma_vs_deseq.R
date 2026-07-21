# ============================================================
#  Limma vs DESeq2 p-value comparison
#  Dataset: WT conditions (FL, HL, ML at T30) vs T0 (control)
# ============================================================


# 0. Load packages ────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("DESeq2", "ggplot2", "ggrepel", "readxl", "dplyr", "patchwork", "ashr", "limma", "edgeR", "tibble", "conflicted")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("DESeq2", "limma", "edgeR")) BiocManager::install(p, ask = FALSE)
    else install.packages(p)
  }
}

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ashr)
library(dplyr)
library(patchwork)

library(limma)
library(edgeR)
library(tibble)


library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


# 1. Load data

data_path <- "Raw_Normalized_Counts_Genes_DESeq2_ToUse.xlsx"   # adjust if needed

raw <- read_excel(data_path, sheet = "All.Genes")


## Round to integers

count_cols <- grep("_Raw.Read.Count$", names(raw), value = TRUE)
counts_mat <- raw[, count_cols] %>%
  mutate(across(everything(), ~ as.integer(round(.)))) %>%
  as.data.frame()

rownames(counts_mat) <- raw$Gene_ID


## Build sample metadata

sample_names <- count_cols
condition <- sub("WT-(FL|HL|ML|T0)-T?30?-[123]_Raw.Read.Count", "\\1", sample_names)
condition[grepl("T0", sample_names)] <- "T0"

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(condition, levels = c("T0", "FL", "HL", "ML")),
  row.names = sample_names
)


# 2 DESeq2 pipeline

## Create DESeq2 object

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ condition
)


## Manual filtering

keepD <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keepD, ]


## Run analysis

dds <- DESeq(dds)


## Extract results for all three contrasts

contrastsD <- list(
  FL_vs_HL = c("condition", "FL", "HL"),
  FL_vs_ML = c("condition", "FL", "ML")
)

results_listD <- lapply(names(contrastsD), function(nm) {
  res <- results(dds,
                 contrastD  = contrastsD[[nm]],
                 alpha     = 0.05,
                 pAdjustMethod = "BH") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene_ID") %>%
    left_join(raw[, c("Gene_ID", "Gene_Name")], by = "Gene_ID") %>%
    mutate(comparison = nm)
  res
})
names(results_listD) <- names(contrastsD)


# 3. Limma pipeline

## Create DGEList

dge <- DGEList(counts = counts_mat, group = col_data$condition)


## Manual filtering

keepL <- rowSums(dge$counts >= 10) >= 3
dge <- dge[keepL, , keepL.lib.sizes = FALSE]


## Normalization

dge <- calcNormFactors(dge, method = "TMM")


## Design matrix, voom, and model fit

design <- model.matrix(~ 0 + condition, data = col_data)
colnames(design) <- levels(col_data$condition)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)


## Contrasts: FL vs HL, FL vs ML

contrast_matrix <- makeContrasts(
  FL_vs_HL = FL - HL,
  FL_vs_ML = FL - ML,
  levels   = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

## Extract results

contrast_names <- colnames(contrast_matrix)

results_listL <- lapply(contrast_names, function(nm) {
  topTable(fit2,
           coef          = nm,
           number        = Inf,
           adjust.method = "BH",
           sort.by       = "P") %>%
    rownames_to_column("Gene_ID") %>%
    rename(
      log2FoldChange = logFC,
      pvalue         = P.Value,
      padj           = adj.P.Val
    ) %>%
    left_join(raw[, c("Gene_ID", "Gene_Name")], by = "Gene_ID") %>%
    mutate(comparison = nm)
})
names(results_listL) <- contrast_names
