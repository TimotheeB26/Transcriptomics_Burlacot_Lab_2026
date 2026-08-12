# Limma-voom differential expression analysis + Venn/Euler diagrams
# Contrasts: DarkvsT0, HLvsT0, FLvsT0

library(readxl)
library(edgeR)
library(limma)
library(eulerr)
library(ggVennDiagram)
library(RColorBrewer)
library(dplyr)
library(conflicted)

conflicts_prefer(base::union)
conflicts_prefer(base::intersect)

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
  HLvsT0   = HL_T30   - T0,
  FLvsT0   = FL_T30   - T0,
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

res_DarkvsT0 <- extract_results(fit2, "DarkvsT0", raw)
res_HLvsT0   <- extract_results(fit2, "HLvsT0",   raw)
res_FLvsT0   <- extract_results(fit2, "FLvsT0",   raw)

# 11. DEG thresholds and gene lists
LFC_THRESH  <- 1
PADJ_THRESH <- 0.05

deg_list <- list(
  `Dark vs T0` = res_DarkvsT0 %>% filter(!is.na(adj.P.Val), adj.P.Val < PADJ_THRESH,
                                          abs(logFC) >= LFC_THRESH) %>% pull(Gene_ID),
  `HL vs T0`   = res_HLvsT0   %>% filter(!is.na(adj.P.Val), adj.P.Val < PADJ_THRESH,
                                          abs(logFC) >= LFC_THRESH) %>% pull(Gene_ID),
  `FL vs T0`   = res_FLvsT0   %>% filter(!is.na(adj.P.Val), adj.P.Val < PADJ_THRESH,
                                          abs(logFC) >= LFC_THRESH) %>% pull(Gene_ID)
)

cat("\nDEG counts per contrast:\n")
print(sapply(deg_list, length))

# 12. Proportional Euler diagram

euler_fit <- euler(deg_list, shape = "ellipse")
venn_cols <- brewer.pal(3, "Set2")

png("limma_Venn_DEGs_proportional.png", width = 2400, height = 2000, res = 300)
plot(euler_fit,
     fills      = list(fill = venn_cols, alpha = 0.55),
     edges      = list(col = venn_cols, lwd = 2),
     labels     = list(font = 2, cex = 1.2, labels = names(deg_list)),
     quantities = list(font = 1, cex = 1.0, type = "counts"),
     main       = list(label = "Venn diagram - limma-voom (|LFC| >= 1, padj < 0.05)",
                        font = 2, cex = 1.3))
dev.off()

# 13. Standard Venn diagram (ggVennDiagram)

p_venn <- ggVennDiagram(
  deg_list,
  label_alpha = 0,
  label       = "count",
  set_color   = venn_cols
) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title    = "limma-voom DEGs: Dark, HL, FL vs T0",
       subtitle = sprintf("|logFC| >= %.0f, FDR < %.2f", LFC_THRESH, PADJ_THRESH),
       fill     = "Gene count") +
  theme_classic() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("limma_Venn_DEGs_ggVennDiagram.png", p_venn, width = 8, height = 7, dpi = 300)

# 14. Intersection summary (7 regions)
cat("\n Intersection summary \n")

Dark <- deg_list[["Dark vs T0"]]
HL   <- deg_list[["HL vs T0"]]
FL   <- deg_list[["FL vs T0"]]

summary_df <- data.frame(
  Region     = c("Dark only", "HL only", "FL only",
                 "Dark Inter HL", "Dark Inter FL", "HL Inter FL",
                 "Dark Inter HL Inter FL"),
  Gene_count = c(
    length(setdiff(Dark, union(HL, FL))),
    length(setdiff(HL, union(Dark, FL))),
    length(setdiff(FL, union(Dark, HL))),
    length(setdiff(intersect(Dark, HL), FL)),
    length(setdiff(intersect(Dark, FL), HL)),
    length(setdiff(intersect(HL, FL), Dark)),
    length(Reduce(intersect, list(Dark, HL, FL)))
  )
)

print(summary_df)
write.csv(summary_df, "limma_DEG_intersection_summary.csv", row.names = FALSE)

# 15. Per-region gene lists, annotated with Gene_Name
dir.create("limma_results", showWarnings = FALSE)

gene_annotation <- raw[, c("Gene_ID", "Gene_Name")]

intersect_genes <- list(
  Dark_only  = setdiff(Dark, union(HL, FL)),
  HL_only    = setdiff(HL, union(Dark, FL)),
  FL_only    = setdiff(FL, union(Dark, HL)),
  Dark_HL    = setdiff(intersect(Dark, HL), FL),
  Dark_FL    = setdiff(intersect(Dark, FL), HL),
  HL_FL      = setdiff(intersect(HL, FL), Dark),
  Dark_HL_FL = Reduce(intersect, list(Dark, HL, FL))
)

for (nm in names(intersect_genes)) {
  ids <- intersect_genes[[nm]]
  if (length(ids) > 0) {
    ann <- gene_annotation %>% filter(Gene_ID %in% ids)
    write.csv(ann, file = sprintf("limma_results/Genes_%s.csv", nm), row.names = FALSE)
  }
}

cat("\n Complete \n")