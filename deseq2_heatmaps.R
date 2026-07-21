# =============================================================================
# DESeq2 Pipeline — WT samples, condition = Tissue + Timepoint
# Groups: FL-T30 | HL-T30 | ML-T30 | T0  (reference = T0)
# Fixes v4:
#   - Significant gene selection now driven by FL-T30 vs HL-T30 contrast
#     (previously: union of FL/HL/ML vs T0)
#   - T0-referenced contrasts (FL/HL/ML vs T0) are still computed and kept
#     for the heatmap display and the final export table, but no longer
#     drive which genes are considered significant
# Fixes v3:
#   - Sample labels: strip Raw.Read.Count from axis titles
#   - HM1: blue-to-red colour scale
#   - HM2: condensed (fixed page size, rasterised, no row names)
#   - HM3: 16 clusters cut from HM2 dendrogram (no re-clustering)
#   - Plot4: faceted horizontal bar for all 16 clusters
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in c("DESeq2","ComplexHeatmap","circlize"))
  if (!requireNamespace(p,quietly=TRUE)) BiocManager::install(p)
for (p in c("readxl","tidyverse","RColorBrewer","cluster","ggplot2","scales"))
  if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(DESeq2); library(ComplexHeatmap); library(circlize)
library(readxl);  library(tidyverse);    library(RColorBrewer)
library(cluster); library(ggplot2);     library(scales)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
FILE_PATH  <- "Raw_Normalized_Counts_Genes_DESeq2_HTSEQ_GO_bis.xlsx"
raw_df     <- read_excel(FILE_PATH, sheet=1)
gene_col   <- colnames(raw_df)[sapply(raw_df, is.character)][1]
go_cols    <- grep("go|term|ontol|annot|function", colnames(raw_df),
                   ignore.case=TRUE, value=TRUE)

# ── 2. COLUMNS ────────────────────────────────────────────────────────────────
count_cols <- c(
  "WT-FL-T30-1_Raw.Read.Count", "WT-FL-T30-2_Raw.Read.Count", "WT-FL-T30-3_Raw.Read.Count",
  "WT-HL-T30-1_Raw.Read.Count", "WT-HL-T30-2_Raw.Read.Count", "WT-HL-T30-3_Raw.Read.Count",
  "WT-ML-T30-1_Raw.Read.Count", "WT-ML-T30-2_Raw.Read.Count", "WT-ML-T30-3_Raw.Read.Count",
  "WT-T0-1_Raw.Read.Count",    "WT-T0-2_Raw.Read.Count",    "WT-T0-3_Raw.Read.Count"
)
# Clean labels shown on heatmap axes (strip "_Raw.Read.Count")
clean_labels <- sub("_Raw\\.Read\\.Count$", "", count_cols)

missing <- setdiff(count_cols, colnames(raw_df))
if (length(missing)>0) stop("Missing columns: ", paste(missing,collapse=", "))

# ── 3. COUNT MATRIX ───────────────────────────────────────────────────────────
count_mat <- raw_df |>
  dplyr::select(all_of(c(gene_col, count_cols))) |>
  column_to_rownames(gene_col) |>
  as.matrix() |> round()
colnames(count_mat) <- clean_labels   # apply clean labels immediately
count_mat[is.na(count_mat)] <- 0
count_mat <- count_mat[rowSums(count_mat)>0, ]

# ── 4. SAMPLE METADATA ────────────────────────────────────────────────────────
parse_condition <- function(nm) {
  nm <- sub("^WT-","",nm); nm <- sub("-[0-9]+$","",nm); nm
}
col_data <- data.frame(
  row.names = clean_labels,
  sample    = clean_labels,
  condition = factor(parse_condition(clean_labels),
                     levels=c("T0","FL-T30","HL-T30","ML-T30"))
)

# ── 5. DESeq2 ─────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=col_data,
                              design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)

# ── 5a. CONTRAST DRIVING SIGNIFICANT GENE SELECTION (FL-T30 vs HL-T30) ───────
res_FLvsHL <- as.data.frame(
  results(dds, contrast=c("condition","FL-T30","HL-T30"), alpha=0.05)
) |> rownames_to_column("Gene")

sig_genes <- res_FLvsHL$Gene[
  !is.na(res_FLvsHL$padj) & res_FLvsHL$padj<0.05 & abs(res_FLvsHL$log2FoldChange)>1
]
message("Significant genes (FL-T30 vs HL-T30, padj<0.05 & |log2FC|>1): ", length(sig_genes))

# ── 5b. T0-REFERENCED CONTRASTS — kept for heatmap context & export table ───
# These no longer determine sig_genes; they are used downstream (heatmaps
# still display all 4 conditions, and the final CSV reports LFC/padj for
# each of FL/HL/ML vs T0 alongside the FL-vs-HL values used for selection).
res_list <- lapply(
  list("FL-T30_vs_T0"=c("condition","FL-T30","T0"),
       "HL-T30_vs_T0"=c("condition","HL-T30","T0"),
       "ML-T30_vs_T0"=c("condition","ML-T30","T0")),
  function(ct) as.data.frame(results(dds,contrast=ct,alpha=0.05)) |>
    rownames_to_column("Gene")
)

go_annot <- raw_df |>
  dplyr::rename(Gene=all_of(gene_col)) |>
  dplyr::select(Gene, all_of(go_cols))

# ── 6. Z-SCORE MATRIX ─────────────────────────────────────────────────────────
z_mat <- assay(vsd)[sig_genes,,drop=FALSE]
z_mat <- t(scale(t(z_mat)))

# ── 7. COLOUR PALETTES ────────────────────────────────────────────────────────
cond_colors <- c("T0"="#4D4D4D","FL-T30"="#1B9E77",
                 "HL-T30"="#D95F02","ML-T30"="#7570B3")
N_CLUSTERS  <- 4
cluster_colors <- setNames(
  colorRampPalette(c(brewer.pal(12,"Paired"), brewer.pal(8,"Dark2")))(N_CLUSTERS),
  as.character(1:N_CLUSTERS))

col_fun_z    <- colorRamp2(c(-2,0,2), c("#3B4CC0","white","#B40426"))

# HM1: blue-to-red for distance (low distance = similar = blue)
dist_range <- range(as.numeric(dist(t(assay(vsd)))))
col_fun_dist <- colorRamp2(
  c(dist_range[1], mean(dist_range), dist_range[2]),
  c("#3B4CC0", "white", "#B40426"))

top_anno <- HeatmapAnnotation(
  Condition=col_data$condition, col=list(Condition=cond_colors),
  annotation_name_side="left")

# ── 8. HEATMAP 1 – SAMPLE DISTANCE (blue → white → red) ──────────────────────
dist_mat <- as.matrix(dist(t(assay(vsd))))

hm1 <- Heatmap(dist_mat, name="Distance", col=col_fun_dist,
               top_annotation=top_anno,
               left_annotation=rowAnnotation(
                 Condition=col_data$condition, col=list(Condition=cond_colors),
                 show_annotation_name=FALSE),
               clustering_distance_rows    = function(m) as.dist(m),
               clustering_distance_columns = function(m) as.dist(m),
               clustering_method_rows="ward.D2", clustering_method_columns="ward.D2",
               show_row_names=TRUE, show_column_names=TRUE,
               row_names_gp=gpar(fontsize=8), column_names_gp=gpar(fontsize=8),
               column_names_rot=45,
               column_title="Sample-to-Sample Distance (VST)",
               column_title_gp=gpar(fontsize=13,fontface="bold"),
               heatmap_legend_param=list(title="Euclidean\nDistance"))

pdf("Heatmap1_Sample_Distance.pdf", width=8, height=7)
draw(hm1, merge_legend=TRUE); dev.off()
message("Saved: Heatmap1_Sample_Distance.pdf")

# ── 9. HEATMAP 2 – Z-SCORE, CONDENSED, SAVE DENDROGRAM ───────────────────────
# Key settings for condensation:
#   - Fixed A4-landscape page (11 x 8.5 in)
#   - use_raster = TRUE always, high quality
#   - row_names hidden (too many genes)
#   - row_dend_width kept small
#   - height argument inside Heatmap controls gene-axis size

hm2 <- Heatmap(z_mat, name="Z-score", col=col_fun_z,
               top_annotation=HeatmapAnnotation(
                 Condition=col_data$condition, col=list(Condition=cond_colors),
                 annotation_name_side="left"),
               clustering_distance_rows="euclidean", clustering_method_rows="ward.D2",
               clustering_distance_columns="euclidean", clustering_method_columns="ward.D2",
               show_row_names=FALSE,          # never show gene names
               show_column_names=TRUE,
               column_names_gp=gpar(fontsize=8), column_names_rot=45,
               row_dend_width=unit(10,"mm"),   # keep dendrogram slim
               column_dend_height=unit(8,"mm"),
               column_title=paste0("z-score - DEGs FL-T30 vs HL-T30 (padj<0.05) - n = ",nrow(z_mat)," genes"),
               column_title_gp=gpar(fontsize=12,fontface="bold"),
               use_raster=TRUE, raster_quality=10,
               height=unit(14,"cm"))           # fixed gene-axis height → fits any n

pdf("Heatmap2_Zscore_AllSigGenes.pdf", width=9, height=8)  # fixed page
hm2_drawn <- draw(hm2, merge_legend=TRUE)
dev.off()
message("Saved: Heatmap2_Zscore_AllSigGenes.pdf")

# ── 10. EXTRACT ROW DENDROGRAM → CUT INTO 16 CLUSTERS ────────────────────────
# We reuse the EXACT same dendrogram that was drawn in HM2 so HM3 clusters
# are directly traceable to the HM2 tree — no second clustering step.

row_dend   <- row_dend(hm2_drawn)           # dendrogram object from drawn hm2
cluster_vec <- cutree(as.hclust(row_dend), k=N_CLUSTERS)
message("Cluster sizes (k=16):")
print(sort(table(cluster_vec), decreasing=TRUE))

# ── 11. HEATMAP 3 – Z-SCORE BY 16 CLUSTERS (from HM2 tree) ───────────────────
# Order rows: keep within-cluster order from original dendrogram
dend_order  <- order.dendrogram(row_dend)    # original leaf order from HM2
cl_ordered  <- cluster_vec[dend_order]        # cluster label in dendrogram order
z_ord       <- z_mat[dend_order, ]
cl_labels   <- as.character(cl_ordered)

hm3 <- Heatmap(z_ord, name="Z-score", col=col_fun_z,
               top_annotation=HeatmapAnnotation(
                 Condition=col_data$condition, col=list(Condition=cond_colors),
                 annotation_name_side="left"),
               left_annotation=rowAnnotation(
                 Cluster=cl_labels, col=list(Cluster=cluster_colors),
                 annotation_name_side="top",
                 annotation_width=unit(3,"mm")),
               cluster_rows=FALSE,               # already ordered from HM2 dendrogram
               clustering_distance_columns="euclidean", clustering_method_columns="ward.D2",
               show_row_names=FALSE,
               show_column_names=TRUE,
               column_names_gp=gpar(fontsize=8), column_names_rot=45,
               column_dend_height=unit(8,"mm"),
               row_split=cl_labels,
               row_title="C%s",                  # short title to save space
               row_title_gp=gpar(fontsize=7,fontface="bold"),
               row_title_rot=0,
               row_gap=unit(0.8,"mm"),           # thin gap between clusters
               column_title=paste0("z-score - DEGs FL-T30 vs HL-T30 - ", N_CLUSTERS, " clusters - n = ",nrow(z_ord)," genes"),
               column_title_gp=gpar(fontsize=12,fontface="bold"),
               use_raster=TRUE, raster_quality=10,
               height=unit(16,"cm"))             # fixed height, fits A4

pdf("Heatmap3_Zscore_9Clusters.pdf", width=10, height=10)  # fixed page
draw(hm3, merge_legend=TRUE); dev.off()
message("Saved: Heatmap3_Zscore_9Clusters.pdf")

# ── 12. PLOT 4 – TOP 8 GO TERMS PAR CLUSTER (POURCENTAGE) ────────────────────
#
# Strategy:
#   A. Parser les GO IDs robustement (séparateurs espace / pipe / ; / ,)
#      et dédupliquer par gène (un gène ne compte qu'une fois par terme)
#   B. Résoudre les noms de termes + ontologie (BP/MF/CC) via GO.db
#   C. Filtrer à une profondeur hiérarchique comparable (depth 3–6) pour
#      comparer des termes au même niveau dans la hiérarchie GO
#   D. Pour chaque cluster : calculer le % de gènes du cluster associés
#      à chaque terme GO, garder le Top 8
#   E. Un PDF par cluster : barplot horizontal coloré par ontologie
# =============================================================================

# ── 12a. Packages ─────────────────────────────────────────────────────────────
for (p in c("GO.db", "AnnotationDbi"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
library(GO.db); library(AnnotationDbi)

# ── 12b. GO IDs → gènes (tidy, dédupliqué) ────────────────────────────────────
go_tidy_raw <- go_annot |>
  dplyr::rename(GO_raw = all_of(go_cols[1])) |>
  filter(!is.na(GO_raw), GO_raw != "") |>
  mutate(GO_id = strsplit(GO_raw, "[\\s|;,]+", perl = TRUE)) |>
  tidyr::unnest(GO_id) |>
  filter(grepl("^GO:\\d{7}$", GO_id)) |>   # GO IDs valides uniquement
  dplyr::select(Gene, GO_id) |>
  distinct()                                # une ligne = (Gene, GO_id) unique

message("Paires GO–gène uniques : ", nrow(go_tidy_raw))

# ── 12c. Noms et ontologie via GO.db ──────────────────────────────────────────
valid_ids <- unique(go_tidy_raw$GO_id)
go_terms  <- AnnotationDbi::select(
  GO.db,
  keys    = valid_ids,
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
) |> dplyr::rename(GO_id = GOID, GO_name = TERM, Ontology = ONTOLOGY) |>
  filter(!is.na(GO_name))

# ── 12d. Profondeur GO (nb d'ancêtres = distance à la racine) ─────────────────
ancestor_env <- list(
  BP = as.list(GOBPANCESTOR),
  MF = as.list(GOMFANCESTOR),
  CC = as.list(GOCCANCESTOR)
)
get_depth <- function(go_id, ont) {
  anc <- ancestor_env[[ont]][[go_id]]
  if (is.null(anc)) return(NA_integer_)
  length(setdiff(anc, "all"))
}
go_terms <- go_terms |>
  rowwise() |>
  mutate(depth = get_depth(GO_id, Ontology)) |>
  ungroup()

# Garder depth 3–6 : termes intermédiaires, comparables entre clusters
# (0–2 = racines trop génériques ; 7+ = feuilles trop spécifiques)
DEPTH_MIN <- 1
DEPTH_MAX <- 9

go_filtered <- go_tidy_raw |>
  inner_join(go_terms, by = "GO_id") |>
  filter(!is.na(depth), depth >= DEPTH_MIN, depth <= DEPTH_MAX)

message(sprintf("Termes GO conservés (depth %d–%d) : %d termes uniques",
                DEPTH_MIN, DEPTH_MAX, n_distinct(go_filtered$GO_id)))

# ── 12e. Associer les clusters aux gènes ──────────────────────────────────────
cluster_df <- data.frame(
  Gene    = names(cluster_vec),
  Cluster = as.character(cluster_vec),
  stringsAsFactors = FALSE
)

go_with_cluster <- go_filtered |>
  inner_join(cluster_df, by = "Gene") |>
  filter(Gene %in% sig_genes)   # restreindre aux gènes DEG significatifs

# ── 12f. Calcul du pourcentage par cluster × terme GO ─────────────────────────
TOP_N <- 8

go_summary <- go_with_cluster |>
  group_by(Cluster, GO_id, GO_name, Ontology) |>
  summarise(n_genes = n_distinct(Gene), .groups = "drop") |>
  mutate(
    GO_label = str_trunc(GO_name, width = 55, side = "right")
  ) |>
  group_by(Cluster) |>
  # On prend le Top N basé sur le nombre de gènes
  slice_max(n_genes, n = TOP_N, with_ties = FALSE) |> 
  ungroup()


# ── 12g. Palette par ontologie ────────────────────────────────────────────────
ont_colors <- c(BP = "#E41A1C", MF = "#377EB8", CC = "#4DAF4A")

# ── 12h. Un PDF par cluster ───────────────────────────────────────────────────
clusters_all <- sort(unique(go_summary$Cluster))

for (cl in clusters_all) {
  
  df_cl <- go_summary |>
    filter(Cluster == cl) |>
    arrange(n_genes) |> # Tri par nombre de gènes
    mutate(GO_label = factor(GO_label, levels = unique(GO_label)))
  
  n_cl_total <- sum(cluster_df$Cluster == cl)
  
  n_cl_total_an <- go_with_cluster |>
    filter(Cluster == cl) |>
    distinct(Gene) |>
    nrow()
  
  if (nrow(df_cl) == 0) {
    message(sprintf("Cluster %s : aucun terme GO — ignoré", cl))
    next
  }
  
  p_go <- ggplot(df_cl, 
                 aes(x = GO_label, y = n_genes, fill = Ontology)) +
    geom_bar(stat = "identity", width = 0.75) +
    # Étiquette affichant uniquement le nombre de gènes
    geom_text(aes(label = n_genes), 
              hjust = -0.3, size = 3, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = ont_colors, 
                      name   = "Ontologie",
                      labels = c(BP = "Biological Process", 
                                 MF = "Molecular Function", 
                                 CC = "Cellular Component")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Espace pour le texte
    coord_flip() +
    labs(
      title    = paste0("Top ", TOP_N, " GO Terms - Cluster ", cl),
      subtitle = sprintf("%d genes (%d annotated)", n_cl_total, n_cl_total_an),
      x = NULL,
      y = "Number of Genes"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      axis.text.y   = element_text(size = 8),
      legend.position = "right"
    )
  
  # Nom de fichier mis à jour
  filename <- paste0("Plot4_GO_Cluster_", cl, "_Top8_Count.pdf")
  h <- max(4, 1.8 + nrow(df_cl) * 0.38)
  pdf(filename, width = 8, height = h)
  print(p_go)
  dev.off()
  message("Saved: ", filename)
}


# ── 12i. Tableau récapitulatif ────────────────────────────────────────────────
write.csv(
  go_summary |>
    # On trie maintenant par Cluster et par nombre de gènes décroissant
    arrange(Cluster, desc(n_genes)) |>
    # On ne sélectionne que les colonnes existantes
    dplyr::select(Cluster, GO_id, GO_name, Ontology, n_genes),
  "Plot4_GO_Top8_Summary.csv",
  row.names = FALSE
)

message("Saved: Plot4_GO_Top8_Summary.csv")


# ── 13. EXPORT TABLE ──────────────────────────────────────────────────────────
res_wide <- Reduce(
  function(a,b) full_join(a,b,by="Gene"),
  lapply(names(res_list), function(nm)
    res_list[[nm]] |>
      dplyr::select(Gene, log2FoldChange, padj) |>
      dplyr::rename_with(function(x) paste0(x,".",nm), -Gene))
)

# Add the FL-T30 vs HL-T30 contrast (the one used to call sig_genes) so the
# selection rationale stays visible alongside the T0-referenced contrasts.
res_wide <- res_wide |>
  full_join(
    res_FLvsHL |>
      dplyr::select(Gene, log2FoldChange, padj) |>
      dplyr::rename_with(function(x) paste0(x,".FL-T30_vs_HL-T30"), -Gene),
    by="Gene"
  )

res_wide |>
  filter(Gene %in% sig_genes) |>
  left_join(data.frame(Gene=names(cluster_vec),
                       Cluster=paste0("C",cluster_vec)), by="Gene") |>
  left_join(go_annot, by="Gene") |>
  arrange(Cluster) |>
  write.csv("DESeq2_sig_genes_clusters.csv", row.names=FALSE)
message("Saved: DESeq2_sig_genes_clusters.csv")

message("\n=== All done! ===")
sessionInfo()