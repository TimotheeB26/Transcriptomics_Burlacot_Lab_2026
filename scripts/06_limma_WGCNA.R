# WGCNA - Weighted Gene Co-expression Network Analysis
# Input:  v$E from 00_limma_voom_normalization.R (voom log2-CPM matrix)
# Design: 6 conditions (T0, Dark_T30, LL_T30, ML_T30, HL_T30, FL_T30) x 3 reps

library(readxl)
library(edgeR)
library(limma)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1. Recreate v$E from the existing limma-voom pipeline
# We re-run the exact same steps as 00_limma_voom_normalization.R rather than
# re-deriving v$E by hand, so filtering/normalization stay identical across
# scripts. If you already have `v` in your session, just skip to section 2.

input_file <- "NRaw_Counts_GO_Terms.xlsx"
raw <- as.data.frame(read_excel(input_file, sheet = "Dataset"))

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

sample_info <- data.frame(sample = colnames(counts), stringsAsFactors = FALSE)
sample_info$group <- sub("_[0-9]+$", "", sample_info$sample)
sample_info$group <- factor(
  sample_info$group,
  levels = c("T0", "Dark_T30", "LL_T30", "ML_T30", "HL_T30", "FL_T30")
)

dge <- DGEList(counts = counts, group = sample_info$group)
keep <- filterByExpr(dge, group = sample_info$group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

v <- voom(dge, design, plot = FALSE)

cat(sprintf("Genes entering WGCNA: %d | Samples: %d\n", nrow(v$E), ncol(v$E)))

# 2. Expression matrix in WGCNA orientation (samples x genes)
datExpr0 <- as.data.frame(t(v$E))

# 3. QC: flag genes/samples with excessive missingness or ~zero variance
# Standard first gate before doing anything else - WGCNA's own recommended
# sanity check. With voom output this should pass cleanly, but it's cheap
# insurance and catches genes that survived filterByExpr with near-zero
# variance across all 6 conditions.
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    cat("Removing", sum(!gsg$goodGenes), "genes failing QC\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("Removing", sum(!gsg$goodSamples), "samples failing QC\n")
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 4. Outlier sample detection via hierarchical clustering
# With only 3 reps/group, a single mislabeled or degraded sample can distort
# module structure disproportionately, so this is worth a visual check rather
# than trusting filterByExpr/voom to have caught it.
sampleTree <- hclust(dist(datExpr0), method = "average")

pdf("WGCNA_01_sample_clustering.pdf", width = 10, height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1, cex.axis = 0.9, cex.main = 1.2)
dev.off()
# No automatic height cut applied here - with n=18 and biologically expected
# condition-driven separation, an automatic cutreeStatic cut risks discarding
# a legitimate condition cluster. Inspect the PDF and set cutHeight manually
# below only if a sample is a clear singleton outlier (uncomment to use):
# clust <- cutreeStatic(sampleTree, cutHeight = <value>, minSize = 10)
# datExpr0 <- datExpr0[clust == 1, ]

datExpr  <- datExpr0
nGenes   <- ncol(datExpr)
nSamples <- nrow(datExpr)

# 5. Trait matrix for module-trait correlation
# Binary indicator ("dummy") encoding of the 6 conditions, since condition is
# categorical/unordered here (T0 is a baseline timepoint, not a light-intensity
# midpoint between Dark and FL). This lets us correlate each module eigengene
# against each condition independently instead of imposing a false ordinal
# light-intensity trait.
sample_info_matched <- sample_info[match(rownames(datExpr), sample_info$sample), ]
traitData <- model.matrix(~ 0 + group, data = sample_info_matched)
colnames(traitData) <- levels(sample_info_matched$group)
rownames(traitData) <- rownames(datExpr)

# 6. Soft-thresholding power selection
# corType = "bicor" (biweight midcorrelation): more robust to outlier samples
# than Pearson, recommended WGCNA default for RNA-seq-derived expression data.
# networkType = "signed": preserves the sign of co-expression (up/up or
# down/down together), which is more biologically interpretable for
# transcriptional modules than "unsigned" (which merges anti-correlated genes
# into the same module) - appropriate here since we care about direction of
# co-regulation across light/dark treatments, not just magnitude of coupling.
powers <- c(1:20)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc      = "bicor",
  verbose     = 5
)

pdf("WGCNA_02_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale free topology R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "blue", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
dev.off()

# Pick the lowest power reaching R^2 >= 0.80 (standard WGCNA threshold for
# signed networks, which need higher powers than unsigned to reach the same
# fit). Falls back to 14 (a common signed-network default for small-to-mid n)
# if no power in the tested range reaches it - inspect the PDF either way
# before trusting this automatically.
sft_row <- sft$fitIndices[sft$fitIndices$SFT.R.sq >= 0.80, ][1, ]
softPower <- if (!is.na(sft_row$Power)) sft_row$Power else 14
cat("Selected soft-thresholding power:", softPower, "\n")

# 7. Network construction + module detection (one-step, blockwise)
# blockwiseModules handles TOM computation, hierarchical clustering, and
# dynamic tree cutting in one call. A single 17k-gene block requires a
# contiguous ~17000x17000 TOM (plus adjacency + working copies) in RAM,
# which commonly fails even on machines with enough *total* free memory -
# so we let WGCNA split genes into blocks (correlation-based preclustering)
# rather than forcing maxBlockSize above nGenes. 5000 is a safe default for
# most machines; drop to 3000-4000 if this still errors, or raise toward
# 8000-10000 if you have >=16GB free and want fewer block-boundary artifacts.
# minModuleSize = 30 and mergeCutHeight = 0.25 are the standard WGCNA starting
# defaults (Langfelder & Horvath); 0.25 merge height means eigengenes
# correlating > 0.75 get merged, which curbs over-fragmentation of modules
# given the relatively small n=18.
net <- blockwiseModules(
  datExpr,
  power              = softPower,
  networkType        = "signed",
  corType            = "bicor",
  TOMType            = "signed",
  maxBlockSize       = 5000,
  minModuleSize      = 30,
  mergeCutHeight     = 0.25,
  numericLabels      = TRUE,
  pamRespectsDendro  = FALSE,
  saveTOMs           = TRUE,
  saveTOMFileBase    = "WGCNA_TOM",
  verbose            = 3
)

cat("Number of blocks used:", length(net$dendrograms), "\n")

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs          <- net$MEs

cat("Modules detected (excluding grey/unassigned):\n")
print(table(moduleColors))

# 8. Dendrogram + module color visualization
# blockwiseModules clusters each block separately, so each block has its own
# dendrogram over its own gene subset (net$blockGenes[[b]]) - plot each one
# rather than assuming a single genome-wide tree.
pdf("WGCNA_03_module_dendrogram.pdf", width = 12, height = 6)
for (b in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[b]]
  plotDendroAndColors(
    net$dendrograms[[b]], moduleColors[block_genes],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = paste("Gene clustering with module assignment - Block", b)
  )
}
dev.off()

# 9. Module-trait relationships
MEs0 <- orderMEs(MEs)
moduleTraitCor <- bicor(MEs0, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Heatmap: module eigengenes x conditions, cell = correlation (annotated p)
pdf("WGCNA_04_module_trait_heatmap.pdf", width = 8, height = max(6, 0.35 * ncol(MEs0)))
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8, 3, 2))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = colnames(traitData),
  yLabels       = colnames(MEs0),
  ySymbols      = colnames(MEs0),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.6,
  zlim          = c(-1, 1),
  main          = "Module-condition relationships (bicor, p-value)"
)
dev.off()

# 10. Hub gene identification per module
# Module membership (kME): correlation of each gene's expression to its
# module eigengene. Genes with high |kME| within their assigned module are
# the most representative "hub" genes - the standard WGCNA definition, more
# informative biologically than raw intramodular connectivity alone since it
# ties back directly to the eigengene summary trait.
geneModuleMembership <- as.data.frame(bicor(datExpr, MEs0, use = "p"))
names(geneModuleMembership) <- paste0("MM.", names(MEs0))

hub_genes <- lapply(unique(moduleColors), function(mod) {
  if (mod == "grey") return(NULL)   # grey = unassigned, not a real module
  me_col <- paste0("ME", moduleLabels[moduleColors == mod][1])
  mm_col <- paste0("MM.", me_col)
  if (!mm_col %in% names(geneModuleMembership)) return(NULL)
  idx <- moduleColors == mod
  mm_vals <- geneModuleMembership[idx, mm_col]
  genes <- colnames(datExpr)[idx]
  ord <- order(-abs(mm_vals))
  data.frame(
    Module  = mod,
    Gene_ID = genes[ord],
    kME     = mm_vals[ord]
  )
})
hub_genes_df <- do.call(rbind, hub_genes)

# Top 10 hub genes per module, annotated with Gene_Name from the original file
top_hubs <- do.call(rbind, lapply(split(hub_genes_df, hub_genes_df$Module), head, 10))
top_hubs <- merge(top_hubs, raw[, c("Gene_ID", "Gene_Name")], by = "Gene_ID", all.x = TRUE)
top_hubs <- top_hubs[order(top_hubs$Module, -abs(top_hubs$kME)), ]

# 11. Export core WGCNA outputs
gene_module_assignment <- data.frame(
  Gene_ID = colnames(datExpr),
  Module_Label = moduleLabels,
  Module_Color = moduleColors
)
gene_module_assignment <- merge(gene_module_assignment, raw[, c("Gene_ID", "Gene_Name")],
                                 by = "Gene_ID", all.x = TRUE)

write.csv(gene_module_assignment, "WGCNA_gene_module_assignment.csv", row.names = FALSE)
write.csv(top_hubs, "WGCNA_top10_hub_genes_per_module.csv", row.names = FALSE)
write.csv(cbind(Condition = colnames(traitData), t(moduleTraitCor)),
          "WGCNA_module_trait_correlations.csv", row.names = FALSE)

save(net, MEs0, moduleColors, moduleLabels, datExpr, traitData,
     file = "WGCNA_network_objects.RData")

cat("\nOutputs written:\n")
cat(" - WGCNA_01_sample_clustering.pdf\n")
cat(" - WGCNA_02_soft_threshold.pdf\n")
cat(" - WGCNA_03_module_dendrogram.pdf\n")
cat(" - WGCNA_04_module_trait_heatmap.pdf\n")
cat(" - WGCNA_gene_module_assignment.csv (", nrow(gene_module_assignment), "genes )\n")
cat(" - WGCNA_top10_hub_genes_per_module.csv\n")
cat(" - WGCNA_module_trait_correlations.csv\n")
cat(" - WGCNA_network_objects.RData (net, MEs0, moduleColors, datExpr, traitData)\n")

# 12. Per-module expression heatmaps (genes x samples), for modules of interest
# This section is self-contained: it reloads the saved .RData objects rather
# than depending on anything computed above, so it can be re-run on its own
# without re-running the (slow) network construction. Just set RDATA_PATH.

RDATA_PATH <- "WGCNA_network_objects.RData"
load(RDATA_PATH)   # restores: net, MEs0, moduleColors, datExpr, traitData

# 12a. limma: DEGs for the FLvsHL contrast
# Re-fits the same voom object with a contrast, since the original script
# only fit the model design, not this specific pairwise contrast. Uses the
# `design` and `v` objects already built in section 1 above (same session) -
# if running this section standalone in a fresh session, re-run section 1
# first so `v`/`design` exist.
fit <- lmFit(v, design)
contrast_FLvsHL <- makeContrasts(FL_T30 - HL_T30, levels = design)
fit2 <- contrasts.fit(fit, contrast_FLvsHL)
fit2 <- eBayes(fit2)

FLvsHL_results <- topTable(fit2, number = Inf, sort.by = "none")
FLvsHL_results$Gene_ID <- rownames(FLvsHL_results)

# Significance cutoff: adj.P.Val < 0.05 AND |logFC| >= 1 (2-fold change).
# Adding the LFC threshold on top of the FDR cutoff keeps statistically
# significant-but-negligible-magnitude changes out of the DEG set - common
# practice when the goal is downstream biological interpretation (e.g. these
# heatmaps) rather than an exhaustive significance list. Adjust LFC_CUTOFF
# if you want a stricter/looser fold-change bar.
LFC_CUTOFF <- 1
DEG_FLvsHL <- FLvsHL_results$Gene_ID[
  FLvsHL_results$adj.P.Val < 0.05 & abs(FLvsHL_results$logFC) >= LFC_CUTOFF
]
cat("FLvsHL DEGs (adj.P.Val < 0.05 & |logFC| >=", LFC_CUTOFF, "):",
    length(DEG_FLvsHL), "of", nrow(FLvsHL_results), "genes tested\n")

write.csv(FLvsHL_results, "WGCNA_FLvsHL_limma_results.csv", row.names = FALSE)

# 12b. Per-module expression heatmaps, restricted to FLvsHL DEGs

# Choose which modules to plot. Default: every non-grey module (grey =
# unassigned genes, not a real co-expression group, so it's excluded).
# To plot only specific modules instead, replace the line below with e.g.:
#   modules_to_plot <- c("turquoise", "brown")
modules_to_plot <- setdiff(unique(moduleColors), "grey")

# Order samples by condition (not alphabetically) so the pattern across
# treatments is easy to read left-to-right.
sample_info_h <- data.frame(
  sample = rownames(traitData),
  group  = sub("_[0-9]+$", "", rownames(traitData)),
  stringsAsFactors = FALSE
)
sample_info_h$group <- factor(
  sample_info_h$group,
  levels = c("T0", "Dark_T30", "LL_T30", "ML_T30", "HL_T30", "FL_T30")
)
sample_info_h <- sample_info_h[order(sample_info_h$group), ]

# Same condition palette/annotation style as the k-means script, so figures
# from both analyses read consistently side by side.
cond_colors <- c(
  "T0"       = "#4D4D4D",
  "Dark_T30" = "#E7298A",
  "LL_T30"   = "#66A61E",
  "ML_T30"   = "#7570B3",
  "HL_T30"   = "#D95F02",
  "FL_T30"   = "#1B9E77"
)
col_fun_z <- colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426"))
top_anno <- HeatmapAnnotation(
  Condition = sample_info_h$group, col = list(Condition = cond_colors),
  annotation_name_side = "left"
)

# Plots one heatmap per module into `pdf_path`. `gene_filter` is either NULL
# (no filtering - full module membership) or a character vector of gene IDs
# to intersect each module's genes against (e.g. DEG_FLvsHL).
plot_module_heatmaps <- function(pdf_path, gene_filter = NULL, subtitle_suffix) {
  pdf(pdf_path, width = 9, height = 8)
  skipped <- character(0)
  for (mod in modules_to_plot) {
    genes_in_module <- colnames(datExpr)[moduleColors == mod]
    if (!is.null(gene_filter)) {
      genes_in_module <- intersect(genes_in_module, gene_filter)
    }

    if (length(genes_in_module) < 2) {  # heatmap needs >=2 rows to cluster
      skipped <- c(skipped, mod)
      next
    }

    mat <- t(datExpr)[genes_in_module, sample_info_h$sample, drop = FALSE]
    mat_z <- t(scale(t(mat)))  # per-gene Z-score across samples, same as col_fun_z's -2/0/2 scale

    hm <- Heatmap(
      mat_z, name = "Z-score", col = col_fun_z,
      top_annotation = top_anno,
      cluster_rows = TRUE, cluster_columns = FALSE,
      column_order = sample_info_h$sample,
      show_row_names = nrow(mat_z) <= 50,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 8),
      column_names_rot = 45,
      column_title = paste0("Module: ", mod, " (n = ", nrow(mat_z), " ", subtitle_suffix, ")"),
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      heatmap_legend_param = list(title = "Z-score")
    )
    draw(hm)
  }
  dev.off()
  skipped
}

# PDF A: unfiltered - full WGCNA module membership
skipped_unfiltered <- plot_module_heatmaps(
  "WGCNA_04a_per_module_heatmaps_all_genes.pdf",
  gene_filter = NULL,
  subtitle_suffix = "genes"
)

# PDF B: filtered - limma -> clustering -> filtering -> heatmap
skipped_filtered <- plot_module_heatmaps(
  "WGCNA_04b_per_module_heatmaps_FLvsHL_DEGs.pdf",
  gene_filter = DEG_FLvsHL,
  subtitle_suffix = "FLvsHL DEGs"
)

plotted_unfiltered <- setdiff(modules_to_plot, skipped_unfiltered)
plotted_filtered   <- setdiff(modules_to_plot, skipped_filtered)
cat(" - WGCNA_FLvsHL_limma_results.csv (", length(DEG_FLvsHL), "DEGs, adj.P.Val < 0.05 )\n")
cat(" - WGCNA_04a_per_module_heatmaps_all_genes.pdf (", length(plotted_unfiltered), "modules plotted )\n")
if (length(skipped_unfiltered) > 0) {
  cat("   Skipped (fewer than 2 genes in module):", paste(skipped_unfiltered, collapse = ", "), "\n")
}
cat(" - WGCNA_04b_per_module_heatmaps_FLvsHL_DEGs.pdf (", length(plotted_filtered), "modules plotted )\n")
if (length(skipped_filtered) > 0) {
  cat("   Skipped (fewer than 2 FLvsHL DEGs in module):", paste(skipped_filtered, collapse = ", "), "\n")
}

cat("\n Complete \n")