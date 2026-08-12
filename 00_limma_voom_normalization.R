# Limma-voom differential expression analysis
# Contrasts: FLvsHL (FL_T30 vs HL_T30), FLvsDark (FL_T30 vs Dark_T30), HLvsDark (HL_T30 vs Dark_T30)

library(readxl)
library(edgeR)
library(limma)

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
  HLvsDark = HL_T30 - Dark_T30,
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
res_HLvsDark <- extract_results(fit2, "HLvsDark", raw)

# 11. Export
write.csv(res_FLvsHL,   "limma_FLvsHL_results.csv",   row.names = FALSE)
write.csv(res_FLvsDark, "limma_FLvsDark_results.csv", row.names = FALSE)
write.csv(res_HLvsDark, "limma_HLvsDark_results.csv", row.names = FALSE)

cat("\n Complete \n")