# ============================================================
#  limma-voom Differential Expression + camera() GO Enrichment
#  Dataset: WT conditions (FL, HL, ML at T30) vs T0 (control)
#  Contrast: FL vs HL
# ============================================================

# ── 0. Install / load packages ────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs_bioc <- c("limma", "edgeR", "GO.db", "AnnotationDbi")
pkgs_cran <- c("ggplot2", "ggrepel", "readxl", "dplyr", "tibble", "stringr", "conflicted")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(limma)
library(edgeR)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)
library(tibble)
library(stringr)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ============================================================
#  PART A — limma-voom differential expression
# ============================================================

# ── A1. Load data ──────────────────────────────────────────────
data_path <- "Raw_Normalized_Counts_Genes_DESeq2_HTSEQ_GO_bis.xlsx"

raw <- read_excel(data_path)

count_cols <- grep("_Raw.Read.Count$", names(raw), value = TRUE)
counts_mat <- raw[, count_cols] %>%
  mutate(across(everything(), ~ as.integer(round(.)))) %>%
  as.data.frame()

rownames(counts_mat) <- raw$Gene_ID

# ── A2. Build sample metadata ──────────────────────────────────
sample_names <- count_cols
condition <- sub("WT-(FL|HL|ML|T0)-T?30?-[123]_Raw.Read.Count", "\\1", sample_names)
condition[grepl("T0", sample_names)] <- "T0"

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(condition, levels = c("T0", "FL", "HL", "ML")),
  row.names = sample_names
)

cat("Sample table:\n")
print(col_data)

# ── A3. Create DGEList, filter, and normalise ──────────────────
dge  <- DGEList(counts = counts_mat, group = col_data$condition)

#keep <- rowSums(dge$counts >= 10) >= 3
keep <- filterByExpr(dge, group = condition, min.count = 10, min.total.count = 15)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("\nGenes after filtering: %d\n", nrow(dge)))

dge <- calcNormFactors(dge, method = "TMM")

# ── A4. Design matrix, voom, and model fit ────────────────────
design <- model.matrix(~ 0 + condition, data = col_data)
colnames(design) <- levels(col_data$condition)   # T0, FL, HL, ML

v    <- voom(dge, design, plot = FALSE)
fit  <- lmFit(v, design)

# ── A5. Contrast: FL vs HL ────────────────────────────────────
contrast_matrix <- makeContrasts(
  FL_vs_HL = FL - HL,
  levels   = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ── A6. Extract results ────────────────────────────────────────
res <- topTable(fit2,
                coef          = "FL_vs_HL",
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
  mutate(comparison = "FL_vs_HL")

write.csv(res,
          file      = "limma_results_FL_vs_HL.csv",
          row.names = FALSE)
cat("Saved: limma_results_FL_vs_HL.csv\n")

# ── A7. Summary of DEGs ────────────────────────────────────────
lfc_cut  <- 1
pval_cut <- 0.05

sig <- res %>% filter(!is.na(padj), padj < pval_cut, abs(log2FoldChange) >= lfc_cut)
up  <- sum(sig$log2FoldChange > 0)
dn  <- sum(sig$log2FoldChange < 0)

cat(sprintf("\n── DEGs FL vs HL (|log2FC| >= %.0f, padj < %.2f) ──\n",
            lfc_cut, pval_cut))
cat(sprintf("  UP (higher in FL): %4d\n", up))
cat(sprintf("  DOWN (lower in FL): %4d\n", dn))
cat(sprintf("  Total:              %4d\n", nrow(sig)))

# ── A8. Volcano plot ───────────────────────────────────────────
vol_colours <- c(
  "NS"               = "grey75",
  "FC only"          = "#FDAE61",
  "Sig only"         = "#ABD9E9",
  "Sig & FC"         = "#D73027"
)

res_vol <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    category = case_when(
      padj < pval_cut & abs(log2FoldChange) >= lfc_cut ~ "Sig & FC",
      padj < pval_cut                                  ~ "Sig only",
      abs(log2FoldChange) >= lfc_cut                   ~ "FC only",
      TRUE                                             ~ "NS"
    )
  )

# Top genes to label: most significant among Sig & FC
top_labels <- res_vol %>%
  filter(category == "Sig & FC") %>%
  slice_min(order_by = padj, n = 40) %>%
  mutate(label = ifelse(!is.na(Gene_Name) & Gene_Name != "", Gene_Name, Gene_ID))

vol_counts       <- table(res_vol$category)
vol_legend_labels <- setNames(
  paste0(names(vol_counts), "  (n=", vol_counts, ")"),
  names(vol_counts)
)

volcano_plot <- ggplot(res_vol, aes(x = log2FoldChange, y = neg_log10_padj,
                                    colour = category)) +
  geom_point(data = ~ filter(., category == "NS"),
             size = 0.7, alpha = 0.35, stroke = 0) +
  geom_point(data = ~ filter(., category != "NS"),
             size = 1.4, alpha = 0.75, stroke = 0) +
  geom_hline(yintercept = -log10(pval_cut),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_text_repel(
    data           = top_labels,
    aes(label      = label),
    size           = 2.5,
    colour         = "black",
    max.overlaps   = 30,
    segment.size   = 0.3,
    segment.colour = "grey50",
    box.padding    = 0.4,
    point.padding  = 0.3
  ) +
  scale_colour_manual(values = vol_colours, labels = vol_legend_labels) +
  labs(
    title    = "Volcano plot - FL vs HL (limma-voom)",
    subtitle = sprintf("padj < %.2f  |  |log2FC| >= %.0f  |  %d genes tested",
                       pval_cut, lfc_cut, nrow(res_vol)),
    x        = expression(log[2]~"FC (FL - HL)"),
    y        = expression(-log[10]~"(adjusted p-value)"),
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
  filename = "Volcano_FL_vs_HL_limma.png",
  plot     = volcano_plot,
  width    = 7, height = 6.5, dpi = 300
)
cat("Saved: Volcano_FL_vs_HL_limma.png\n")

cat("\n Part A (limma-voom DE) complete.\n")


# ============================================================
#  PART B — GO Term Enrichment using limma::camera()
#  GO annotations: GO_Terms_Mart column (space-separated GO IDs)
#  GO hierarchy propagated to ancestor terms, per ontology (BP/MF/CC)
# ============================================================

# ── B1. Config ────────────────────────────────────────────────
go_annot_path  <- data_path  # same file as Part A; has GO_Terms_Mart column
minGSSize      <- 5      # min genes per GO set (within tested universe)
maxGSSize      <- 500    # max genes per GO set
# inter.gene.cor: NA -> estimate per set (strict error control, per camera() docs);
#                 0.01 -> preset, favors ranking biologically-interpretable sets higher.
inter_gene_cor <- 0.01
fdr_cutoff     <- 0.05

# v, design, contrast_matrix all come from Part A above — no reload needed.

# ── B2. Load GO annotations & build gene -> direct GO map ─────
go_raw <- read_excel(go_annot_path) %>%
  dplyr::select(Gene_ID, GO_Terms_Mart) %>%
  dplyr::filter(!is.na(GO_Terms_Mart), GO_Terms_Mart != "")

gene2go_direct <- setNames(
  str_split(go_raw$GO_Terms_Mart, "\\s+"),
  go_raw$Gene_ID
)

cat(sprintf("Genes with direct GO annotation: %d\n", length(gene2go_direct)))

# ── B3. Ontology (BP/MF/CC) + term name for every observed GO ID ──
all_go_ids <- unique(unlist(gene2go_direct))

go_onto <- AnnotationDbi::select(
  GO.db,
  keys    = all_go_ids,
  columns = c("ONTOLOGY", "TERM"),
  keytype = "GOID"
) %>%
  dplyr::filter(!is.na(ONTOLOGY))   # drops obsolete / unmapped GO IDs

id2onto <- setNames(go_onto$ONTOLOGY, go_onto$GOID)
id2term <- setNames(go_onto$TERM,     go_onto$GOID)

# ── B4. Propagate GO hierarchy: gene -> direct + all ancestor terms ──
# Expand ancestors ONCE per unique GO term (not per gene), then do fast
# list lookups for each gene. Avoids thousands of redundant mget() calls.
cat("Propagating GO hierarchy (BP/MF/CC ancestors)...\n")

ancestor_dbs <- list(BP = GOBPANCESTOR, MF = GOMFANCESTOR, CC = GOCCANCESTOR)

# Build term -> (self + all ancestors) lookup, one ontology at a time
term2full <- list()
for (onto in names(ancestor_dbs)) {
  terms_onto <- names(id2onto)[id2onto == onto]
  if (length(terms_onto) == 0) next

  db      <- ancestor_dbs[[onto]]
  mapped  <- terms_onto[terms_onto %in% mappedkeys(db)]
  unmapped <- setdiff(terms_onto, mapped)

  if (length(mapped) > 0) {
    anc_list <- AnnotationDbi::mget(mapped, db)                 # one call for all terms in this ontology
    anc_list <- lapply(anc_list, function(x) setdiff(x, c(NA, "all")))
    term2full[mapped] <- Map(function(trm, anc) unique(c(trm, anc)),
                              mapped, anc_list)
  }
  # Terms with no ancestor mapping (e.g. root terms) just map to themselves
  if (length(unmapped) > 0) {
    term2full[unmapped] <- as.list(unmapped)
  }
}

# Fast per-gene lookup: union of term2full[[t]] for each direct term t
propagate_gene_go <- function(direct_ids) {
  direct_ids <- direct_ids[direct_ids %in% names(term2full)]
  if (length(direct_ids) == 0) return(character(0))
  unique(unlist(term2full[direct_ids], use.names = FALSE))
}

gene2go_full <- lapply(gene2go_direct, propagate_gene_go)

# ── B5. Restrict to genes actually tested in the voom/limma model ────
universe     <- rownames(v$E)
gene2go_full <- gene2go_full[names(gene2go_full) %in% universe]

go2genes <- split(
  rep(names(gene2go_full), lengths(gene2go_full)),
  unlist(gene2go_full)
)

# Gene-set size filter (applied to the size WITHIN the tested universe)

cat(sprintf("GO terms retained before size filter (%d-%d genes): %d\n",
            minGSSize, maxGSSize, length(go2genes)))

go_sizes <- lengths(go2genes)
go2genes <- go2genes[go_sizes >= minGSSize & go_sizes <= maxGSSize]
go2genes <- go2genes[names(go2genes) %in% names(id2onto)]  # keep only terms w/ known ontology

cat(sprintf("GO terms retained after size filter (%d-%d genes): %d\n",
            minGSSize, maxGSSize, length(go2genes)))

# camera()'s `index` wants row positions in y, not gene IDs
go2index <- lapply(go2genes, function(g) match(g, universe))

# ── B6. Run camera(), separately per ontology (independent BH correction) ──
ontologies       <- c("BP", "MF", "CC")
camera_res_list  <- list()

for (onto in ontologies) {
  onto_terms <- names(go2index)[id2onto[names(go2index)] == onto]
  if (length(onto_terms) == 0) next

  res_onto <- camera(
    y              = v,
    index          = go2index[onto_terms],
    design         = design,
    contrast       = contrast_matrix[, "FL_vs_HL"],
    inter.gene.cor = inter_gene_cor,
    sort           = FALSE
  ) %>%
    rownames_to_column("GO_ID") %>%
    mutate(
      Ontology = onto,
      Term     = id2term[GO_ID],
      FDR      = p.adjust(PValue, method = "BH")   # per-ontology BH correction
    ) %>%
    arrange(PValue)

  camera_res_list[[onto]] <- res_onto
  cat(sprintf("  %s: %d terms tested\n", onto, nrow(res_onto)))
}

camera_results <- bind_rows(camera_res_list) %>%
  dplyr::select(GO_ID, Term, Ontology, NGenes, any_of("Correlation"),
                Direction, PValue, FDR) %>%
  arrange(Ontology, PValue)

# ── B7. Export ─────────────────────────────────────────────────
write.csv(camera_results, "camera_GO_enrichment_FL_vs_HL.csv", row.names = FALSE)
cat("Saved: camera_GO_enrichment_FL_vs_HL.csv\n")

sig_results <- camera_results %>% dplyr::filter(FDR < fdr_cutoff)
cat(sprintf("\nSignificant GO terms (FDR < %.2f): %d\n", fdr_cutoff, nrow(sig_results)))
print(table(sig_results$Ontology))

# ── B8. Dotplot of top significant terms per ontology (separate pages, same PDF) ──
top_per_onto <- camera_results %>%
  dplyr::filter(FDR < fdr_cutoff) %>%
  group_by(Ontology) %>%
  slice_min(order_by = FDR, n = 30) %>%
  ungroup() %>%
  mutate(
    Term_wrapped = str_wrap(Term, width = 45),
    Term_wrapped = factor(Term_wrapped, levels = rev(unique(Term_wrapped)))
  )

if (nrow(top_per_onto) > 0) {

  make_onto_plot <- function(df, onto_name) {
    ggplot(df, aes(x = -log10(FDR), y = Term_wrapped,
                    colour = Direction, size = NGenes)) +
      geom_point() +
      scale_colour_manual(values = c(Up = "#D73027", Down = "#4575B4")) +
      labs(
        title = sprintf("Top GO terms - %s - camera (FLvsHL)", onto_name),
        x     = expression(-log[10]~"(FDR)"),
        y     = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.y   = element_text(size = 9),
        axis.text.x   = element_text(size = 10),
        legend.text   = element_text(size = 10),
        legend.title  = element_text(size = 11),
        plot.title    = element_text(face = "bold", size = 14),
        plot.margin   = margin(10, 15, 10, 10)
      )
  }

  onto_list <- split(top_per_onto, top_per_onto$Ontology)

  # A single PDF device has one fixed page size, so size it to the
  # ontology with the most terms; smaller-ontology pages just get
  # a bit of blank space below the last point, which is fine.
  page_height <- max(6, sapply(onto_list, function(df) 0.35 * nrow(df) + 1.5))

  pdf_file <- "camera_GO_dotplot_FL_vs_HL.pdf"
  pdf(pdf_file, width = 10, height = page_height, onefile = TRUE)
  for (onto_name in names(onto_list)) {
    print(make_onto_plot(onto_list[[onto_name]], onto_name))
  }
  dev.off()

  cat(sprintf("Saved: %s (%d ontology pages, %.1f x %.1f in)\n",
              pdf_file, length(onto_list), 10, page_height))
} else {
  cat("No significant terms at the chosen FDR cutoff — dotplot skipped.\n")
}

cat("\n Part B (camera GO enrichment) complete.\n")
cat("\nAll done.\n")
