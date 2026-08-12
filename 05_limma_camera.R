# limma-voom differential expression + camera() GO enrichment
# Dataset: WT conditions (Dark, LL, ML, HL, FL at T30) vs T0 (control)
# Contrasts: FL vs HL, FL vs Dark, HL vs Dark, LL vs Dark, ML vs Dark

# 1. Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs_bioc <- c("limma", "edgeR", "GO.db", "AnnotationDbi")
pkgs_cran <- c("readxl", "ggplot2", "ggrepel", "dplyr", "tibble", "stringr", "conflicted")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(edgeR)
library(limma)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(stringr)
library(conflicted)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::setdiff)

# 2. Output directories
out_dir      <- "limma_camera_outputs"
dir_go       <- file.path(out_dir, "camera_GO_enrichment")
dir_violin   <- file.path(out_dir, "LFC_violin_plots")
dir_combined <- file.path(out_dir, "combined_dotplots")

for (d in c(out_dir, dir_go, dir_violin, dir_combined)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 3. Import data
input_file <- "NRaw_Counts_GO_Terms.xlsx"
raw <- as.data.frame(read_excel(input_file, sheet = "Dataset"))

# 4. Build the count matrix
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

# 5. Sample metadata
sample_info <- data.frame(sample = colnames(counts), stringsAsFactors = FALSE)
sample_info$group <- sub("_[0-9]+$", "", sample_info$sample)
sample_info$group <- factor(
  sample_info$group,
  levels = c("T0", "Dark_T30", "LL_T30", "ML_T30", "HL_T30", "FL_T30")
)

# 6. DGEList + low-count filtering
dge <- DGEList(counts = counts, group = sample_info$group)

keep <- filterByExpr(dge, group = sample_info$group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("Genes retained after filterByExpr: %d / %d\n", sum(keep), length(keep)))

# 7. TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# 8. Design matrix
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

# 9. Voom transformation
v <- voom(dge, design, plot = FALSE)

# 10. Linear model fit
fit <- lmFit(v, design)

# 11. Contrasts
contrast_matrix <- makeContrasts(
  FLvsHL   = FL_T30 - HL_T30,
  FLvsDark = FL_T30 - Dark_T30,
  HLvsDark = HL_T30 - Dark_T30,
  LLvsDark = LL_T30 - Dark_T30,
  MLvsDark = ML_T30 - Dark_T30,
  levels   = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 12. Contrast definitions used throughout the per-contrast pipeline and Part C
contrast_defs <- list(
  FLvsHL   = list(coef = "FLvsHL",   label = "FL vs HL",   axis_label = "FL - HL"),
  FLvsDark = list(coef = "FLvsDark", label = "FL vs Dark", axis_label = "FL - Dark"),
  HLvsDark = list(coef = "HLvsDark", label = "HL vs Dark", axis_label = "HL - Dark"),
  LLvsDark = list(coef = "LLvsDark", label = "LL vs Dark", axis_label = "LL - Dark"),
  MLvsDark = list(coef = "MLvsDark", label = "ML vs Dark", axis_label = "ML - Dark")
)

# Contrasts (and their x-axis order) used in the combined multi-contrast
# GO dotplot built in Part C below.
combined_plot_contrasts <- c("LLvsDark", "MLvsDark", "HLvsDark", "FLvsDark")

# 13. Extract per-contrast DE results right after the contrasts are defined
# (same approach as 00_limma_voom_normalization.R's extract_results()).
extract_results <- function(fit_obj, coef_name, tag, annotation) {
  topTable(fit_obj, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "P") %>%
    rownames_to_column("Gene_ID") %>%
    dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val) %>%
    left_join(annotation[, c("Gene_ID", "Gene_Name")], by = "Gene_ID") %>%
    mutate(comparison = tag)
}

res_list <- setNames(
  lapply(names(contrast_defs), function(nm) {
    extract_results(fit2, contrast_defs[[nm]]$coef, nm, raw)
  }),
  names(contrast_defs)
)

# 14. GO enrichment config (shared across contrasts)
go_annot_path  <- input_file  # same file as Part A; has GO_Terms column
minGSSize      <- 5      # min genes per GO set (within tested universe)
maxGSSize      <- 500    # max genes per GO set
fdr_cutoff     <- 0.05
lfc_cut        <- 1
pval_cut       <- 0.05
inter_gene_cor <- 0.01 # 0.01 recommended in camera documentation

terms_of_interest <- c("translation", "signal transduction")   # edit as needed

# 15. Load GO annotations & build gene -> direct GO map (shared)
go_raw <- read_excel(go_annot_path) %>%
  dplyr::select(Gene_ID, GO_Terms) %>%
  dplyr::filter(!is.na(GO_Terms), GO_Terms != "")

gene2go_direct <- setNames(
  str_split(go_raw$GO_Terms, "\\s+"),
  go_raw$Gene_ID
)

cat(sprintf("Genes with direct GO annotation: %d\n", length(gene2go_direct)))

# 16. Ontology (BP/MF/CC) + term name for every observed GO ID (shared)
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

# 17. Propagate GO hierarchy: gene -> direct + all ancestor terms (shared)
cat("Propagating GO hierarchy (BP/MF/CC ancestors)...\n")

ancestor_dbs <- list(BP = GOBPANCESTOR, MF = GOMFANCESTOR, CC = GOCCANCESTOR)

term2full <- list()
for (onto in names(ancestor_dbs)) {
  terms_onto <- names(id2onto)[id2onto == onto]
  if (length(terms_onto) == 0) next

  db      <- ancestor_dbs[[onto]]
  mapped  <- terms_onto[terms_onto %in% mappedkeys(db)]
  unmapped <- setdiff(terms_onto, mapped)

  if (length(mapped) > 0) {
    anc_list <- AnnotationDbi::mget(mapped, db)
    anc_list <- lapply(anc_list, function(x) setdiff(x, c(NA, "all")))
    term2full[mapped] <- Map(function(trm, anc) unique(c(trm, anc)),
                             mapped, anc_list)
  }
  if (length(unmapped) > 0) {
    term2full[unmapped] <- as.list(unmapped)
  }
}

propagate_gene_go <- function(direct_ids) {
  direct_ids <- direct_ids[direct_ids %in% names(term2full)]
  if (length(direct_ids) == 0) return(character(0))
  unique(unlist(term2full[direct_ids], use.names = FALSE))
}

gene2go_full <- lapply(gene2go_direct, propagate_gene_go)

# 18. Restrict to genes actually tested in the voom/limma model (shared)
universe     <- rownames(v$E)
gene2go_full <- gene2go_full[names(gene2go_full) %in% universe]

go2genes_all <- split(
  rep(names(gene2go_full), lengths(gene2go_full)),
  unlist(gene2go_full)
)

go2genes <- go2genes_all

cat(sprintf("GO terms retained before size filter (%d-%d genes): %d\n",
            minGSSize, maxGSSize, length(go2genes)))

go_sizes <- lengths(go2genes)
go2genes <- go2genes[go_sizes >= minGSSize & go_sizes <= maxGSSize]
go2genes <- go2genes[names(go2genes) %in% names(id2onto)]  # keep only terms w/ known ontology

cat(sprintf("GO terms retained after size filter (%d-%d genes): %d\n",
            minGSSize, maxGSSize, length(go2genes)))

go2index <- lapply(go2genes, function(g) match(g, universe))

ontologies <- c("BP", "MF", "CC")

# 19. Helper: most-general (parent-only) filter for a set of significant GO ids
filter_most_general_terms <- function(go_ids, ancestor_db) {
  if (length(go_ids) <= 1) return(go_ids)

  mapped   <- go_ids[go_ids %in% mappedkeys(ancestor_db)]
  anc_list <- if (length(mapped) > 0) AnnotationDbi::mget(mapped, ancestor_db) else list()

  is_child_of_another <- function(term) {
    ancestors_of_term <- anc_list[[term]]
    if (is.null(ancestors_of_term)) return(FALSE)
    any(go_ids[go_ids != term] %in% ancestors_of_term)
  }

  keep <- !vapply(go_ids, is_child_of_another, logical(1))
  go_ids[keep]
}

# 20. Helper: dotplot builder for top significant terms per ontology
make_onto_plot <- function(df, onto_name, contrast_label) {
  ggplot(df, aes(x = -log10(FDR), y = Term_wrapped,
                 colour = Direction, size = NGenes)) +
    geom_point() +
    scale_colour_manual(values = c(Up = "#D73027", Down = "#4575B4")) +
    labs(
      title = sprintf("Top GO terms - %s - camera (%s)", onto_name, contrast_label),
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

# 21. Helper: combined multi-contrast dotplot (one plot per ontology)
# x = contrast, y = GO term, colour = Direction, size = -log10(FDR).
# Includes every term that is significant (FDR < fdr_cutoff) in at least
# one of `contrast_order`; each such term is plotted across ALL of those
# contrasts (even where it isn't itself significant there) so the pattern
# across conditions is visible.
make_combined_dotplot <- function(all_camera_results, onto_name, contrast_order, contrast_labels,
                                   go_id_filter = NULL, title_label = NULL) {

  df_onto <- all_camera_results %>%
    dplyr::filter(Ontology == onto_name, Contrast %in% contrast_order)

  # Optional restriction to a specific set of GO IDs (e.g. most-general BP terms)
  if (!is.null(go_id_filter)) {
    df_onto <- df_onto %>% dplyr::filter(GO_ID %in% go_id_filter)
  }

  if (is.null(title_label)) title_label <- onto_name

  sig_terms <- df_onto %>%
    dplyr::filter(FDR < fdr_cutoff) %>%
    dplyr::pull(GO_ID) %>%
    unique()

  if (length(sig_terms) == 0) return(NULL)

  # Keep a dot only where the term is itself enriched (FDR < fdr_cutoff) in
  # that specific contrast — non-significant cells for an otherwise-included
  # term are left blank rather than plotted.
  plot_df <- df_onto %>%
    dplyr::filter(GO_ID %in% sig_terms, FDR < fdr_cutoff) %>%
    mutate(
      Contrast     = factor(Contrast, levels = contrast_order, labels = contrast_labels),
      Term_wrapped = str_wrap(Term, width = 45)
    )

  # Order terms top-to-bottom by their best (minimum) FDR across the shown contrasts
  term_order <- plot_df %>%
    group_by(Term_wrapped) %>%
    summarise(best_fdr = min(FDR), .groups = "drop") %>%
    arrange(best_fdr) %>%
    pull(Term_wrapped)

  plot_df <- plot_df %>%
    mutate(Term_wrapped = factor(Term_wrapped, levels = rev(term_order)))

  n_terms <- length(term_order)

  ggplot(plot_df, aes(x = Contrast, y = Term_wrapped,
                       colour = Direction, size = -log10(FDR))) +
    geom_point() +
    scale_colour_manual(values = c(Up = "#D73027", Down = "#4575B4")) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)")) +
    labs(
      title = sprintf("GO enrichment across contrasts - %s (camera)", title_label),
      x     = NULL,
      y     = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x      = element_text(size = 11, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 9),
      legend.text      = element_text(size = 10),
      legend.title     = element_text(size = 11),
      plot.title       = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(10, 15, 10, 10)
    ) -> p

  list(plot = p, n_terms = n_terms)
}

# 22. Helper: violin/box/jitter plot of log2FC for a set of enriched terms
# The jittered points are drawn AFTER the boxplot so the red dots sit on
# top of the IQR box instead of being buried underneath it.
make_all_terms_lfc_plot <- function(onto_name, sig_onto, go2genes, gene_lfc_lookup,
                                     contrast_label, axis_label) {

  term_gene_df <- lapply(seq_len(nrow(sig_onto)), function(i) {
    go_id <- sig_onto$GO_ID[i]
    genes <- go2genes[[go_id]]
    if (is.null(genes) || length(genes) == 0) return(NULL)
    data.frame(
      GO_ID   = go_id,
      Term    = sig_onto$Term[i],
      FDR     = sig_onto$FDR[i],
      Gene_ID = genes,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    left_join(gene_lfc_lookup, by = "Gene_ID") %>%
    dplyr::filter(!is.na(log2FoldChange))

  if (nrow(term_gene_df) == 0) return(NULL)

  term_order <- sig_onto %>% arrange(FDR) %>% pull(Term) %>% unique()
  term_gene_df <- term_gene_df %>%
    mutate(
      Term_wrapped = str_wrap(Term, width = 30),
      Term_wrapped = factor(Term_wrapped,
                            levels = rev(str_wrap(term_order, width = 30) %>% unique()))
    )

  n_by_term <- term_gene_df %>%
    group_by(Term_wrapped) %>%
    summarise(n_genes = dplyr::n(), max_lfc = max(log2FoldChange), .groups = "drop")

  n_terms <- length(unique(term_gene_df$Term_wrapped))

  p <- ggplot(term_gene_df, aes(x = log2FoldChange, y = Term_wrapped)) +
    geom_violin(
      fill        = "grey85",
      colour      = "grey40",
      width       = 0.8,
      trim        = FALSE,
      scale       = "width"
    ) +
    geom_boxplot(width = 0.18, outlier.shape = NA, colour = "black", fill = "white", linewidth = 0.6) +
    geom_jitter(height = 0.05, width = 0, size = 1.2, alpha = 0.5, colour = "#D73027") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_text(
      data        = n_by_term,
      aes(x = max_lfc, y = Term_wrapped, label = sprintf("n=%d", n_genes)),
      hjust       = -0.2,
      size        = 2.6,
      inherit.aes = FALSE
    ) +
    labs(
      title    = sprintf("Log2FC distribution by enriched GO Term (Limma-Camera) - %s - %s", onto_name, contrast_label),
      subtitle = "Dots: genes - Violins: distribution (inner box: median/IQR)",
      y        = NULL,
      x        = bquote(log[2]~"FC ("*.(axis_label)*")")
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.y      = element_text(size = 8),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(10, 15, 10, 15)
    )

  list(plot = p, n_terms = n_terms)
}

# 23. Per-contrast pipeline: camera GO enrichment, dotplots, and log2FC violin plots
# (DE results themselves were already extracted in step 13)
run_contrast_analysis <- function(contrast_name, cdef) {

  coef_name      <- cdef$coef
  contrast_label <- cdef$label     # e.g. "FL vs HL"
  axis_label     <- cdef$axis_label # e.g. "FL - HL"
  tag            <- contrast_name  # e.g. "FLvsHL" -- used in filenames

  cat(sprintf("\n\n================ Running contrast: %s ================\n", contrast_label))

  # 23a. Retrieve results (already extracted in step 13, right after the contrasts were defined)
  res <- res_list[[contrast_name]]

  # 23b. Summary of DEGs
  sig <- res %>% filter(!is.na(padj), padj < pval_cut, abs(log2FoldChange) >= lfc_cut)
  up  <- sum(sig$log2FoldChange > 0)
  dn  <- sum(sig$log2FoldChange < 0)

  cat(sprintf("\n── DEGs %s (|log2FC| >= %.0f, padj < %.2f) ──\n",
              contrast_label, lfc_cut, pval_cut))
  cat(sprintf("  UP: %4d\n", up))
  cat(sprintf("  DOWN: %4d\n", dn))
  cat(sprintf("  Total:              %4d\n", nrow(sig)))

  cat(sprintf("\n Part A (limma-voom DE) complete for %s.\n", contrast_label))

  # ============================================================
  #  PART B — GO Term Enrichment using limma::camera()
  # ============================================================

  # 23c. Run camera(), separately per ontology (independent BH correction)
  camera_res_list <- list()

  for (onto in ontologies) {
    onto_terms <- names(go2index)[id2onto[names(go2index)] == onto]
    if (length(onto_terms) == 0) next

    res_onto <- camera(
      y              = v,
      index          = go2index[onto_terms],
      design         = design,
      contrast       = contrast_matrix[, coef_name],
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

  # 23d. Export camera results
  camera_csv <- file.path(dir_go, sprintf("camera_GO_enrichment_%s.csv", tag))
  write.csv(camera_results, camera_csv, row.names = FALSE)
  cat(sprintf("Saved: %s\n", camera_csv))

  sig_results <- camera_results %>% dplyr::filter(FDR < fdr_cutoff)
  cat(sprintf("\nSignificant GO terms (FDR < %.2f): %d\n", fdr_cutoff, nrow(sig_results)))
  print(table(sig_results$Ontology))

  # 23e. Dotplot of top significant terms per ontology (separate pages, same PDF)
  top_per_onto <- camera_results %>%
    dplyr::filter(FDR < fdr_cutoff) %>%
    group_by(Ontology) %>%
    slice_min(order_by = FDR, n = 300) %>%
    ungroup() %>%
    mutate(
      Term_wrapped = str_wrap(Term, width = 45),
      Term_wrapped = factor(Term_wrapped, levels = rev(unique(Term_wrapped)))
    )

  if (nrow(top_per_onto) > 0) {

    onto_list <- split(top_per_onto, top_per_onto$Ontology)

    page_height <- max(6, sapply(onto_list, function(df) 0.35 * nrow(df) + 1.5))

    pdf_file <- file.path(dir_go, sprintf("camera_GO_dotplot_%s.pdf", tag))
    pdf(pdf_file, width = 10, height = page_height, onefile = TRUE)
    for (onto_name in names(onto_list)) {
      print(make_onto_plot(onto_list[[onto_name]], onto_name, contrast_label))
    }
    dev.off()

    cat(sprintf("Saved: %s (%d ontology pages, %.1f x %.1f in)\n",
                pdf_file, length(onto_list), 10, page_height))
  } else {
    cat("No significant terms at the chosen FDR cutoff — dotplot skipped.\n")
  }

  # 23f. BP plot keeping only the most general term per hierarchy
  bp_sig_all <- sig_results %>% dplyr::filter(Ontology == "BP")

  if (nrow(bp_sig_all) > 0) {

    general_bp_ids <- filter_most_general_terms(bp_sig_all$GO_ID, GOBPANCESTOR)

    cat(sprintf("BP significant terms: %d total -> %d most-general (parents only, children dropped)\n",
                nrow(bp_sig_all), length(general_bp_ids)))

    bp_general <- bp_sig_all %>%
      dplyr::filter(GO_ID %in% general_bp_ids) %>%
      slice_min(order_by = FDR, n = 300) %>%
      mutate(
        Term_wrapped = str_wrap(Term, width = 45),
        Term_wrapped = factor(Term_wrapped, levels = rev(unique(Term_wrapped)))
      )

    bp_general_plot <- make_onto_plot(bp_general, "BP - most general terms only", contrast_label)

    bp_general_file <- file.path(dir_go, sprintf("camera_GO_dotplot_BP_most_general_%s.png", tag))
    ggsave(
      filename = bp_general_file,
      plot     = bp_general_plot,
      width    = 10,
      height   = max(4, 0.35 * nrow(bp_general) + 1.5),
      dpi      = 300
    )
    cat(sprintf("Saved: %s\n", bp_general_file))

  } else {
    cat("No significant BP terms — most-general BP plot skipped.\n")
  }

  # 23g. Log2FC distribution per GO term (hand-picked terms_of_interest)
  # Jitter drawn after the box so the red dots sit on top of the IQR box.
  term_ids_tbl <- AnnotationDbi::select(
    GO.db,
    keys    = terms_of_interest,
    keytype = "TERM",
    columns = c("GOID", "TERM")
  )

  lfc_by_term <- lapply(terms_of_interest, function(tn) {
    ids   <- term_ids_tbl$GOID[term_ids_tbl$TERM == tn]
    ids   <- ids[!is.na(ids)]
    genes <- unique(unlist(go2genes_all[ids]))
    if (length(genes) == 0) return(NULL)
    data.frame(GO_Term = tn, Gene_ID = genes, stringsAsFactors = FALSE)
  }) %>%
    bind_rows() %>%
    left_join(res[, c("Gene_ID", "log2FoldChange")], by = "Gene_ID") %>%
    dplyr::filter(!is.na(log2FoldChange))

  if (nrow(lfc_by_term) > 0) {

    n_by_term <- lfc_by_term %>%
      group_by(GO_Term) %>%
      summarise(n_genes = dplyr::n(), max_lfc = max(log2FoldChange), .groups = "drop")

    missing_terms <- setdiff(terms_of_interest, unique(lfc_by_term$GO_Term))
    if (length(missing_terms) > 0) {
      cat(sprintf("No annotated/tested genes found for: %s\n", paste(missing_terms, collapse = ", ")))
    }

    lfc_term_plot <- ggplot(lfc_by_term, aes(x = GO_Term, y = log2FoldChange)) +
      geom_violin(
        fill        = "grey85",
        colour      = "grey40",
        width       = 0.7,
        trim        = FALSE,
        scale       = "width"
      ) +
      geom_boxplot(width = 0.15, outlier.shape = NA, colour = "black", fill = "white", linewidth = 0.6) +
      geom_jitter(width = 0.05, height = 0, size = 1.6, alpha = 0.6, colour = "#D73027") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
      geom_text(
        data        = n_by_term,
        aes(x = GO_Term, y = max_lfc, label = sprintf("n=%d", n_genes)),
        vjust       = -0.8,
        size        = 3,
        inherit.aes = FALSE
      ) +
      labs(
        title    = sprintf("Log2FC distribution by GO Term (%s)", contrast_label),
        subtitle = "Violin = distribution of log2FC across genes annotated to the term; dots = individual genes",
        x        = NULL,
        y        = bquote(log[2]~"FC ("*.(axis_label)*")")
      ) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x      = element_text(size = 11),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank()
      )

    lfc_term_file <- file.path(dir_violin, sprintf("GO_term_LFC_violin_%s.png", tag))
    ggsave(
      filename = lfc_term_file,
      plot     = lfc_term_plot,
      width    = max(6, 1.5 * length(unique(lfc_by_term$GO_Term)) + 3),
      height   = 6,
      dpi      = 300
    )
    cat(sprintf("Saved: %s\n", lfc_term_file))

  } else {
    cat("No genes found for the specified GO terms — LFC violin plot skipped.\n")
  }

  # 23h. Log2FC distribution violin + dots for all enriched GO terms, one plot per ontology
  gene_lfc_lookup <- res %>% dplyr::select(Gene_ID, log2FoldChange)

  for (onto in ontologies) {

    sig_onto <- sig_results %>% dplyr::filter(Ontology == onto) %>% arrange(FDR)

    if (nrow(sig_onto) == 0) {
      cat(sprintf("No significant %s terms — all-terms LFC violin plot skipped.\n", onto))
      next
    }

    result <- make_all_terms_lfc_plot(onto, sig_onto, go2genes, gene_lfc_lookup,
                                       contrast_label, axis_label)

    if (is.null(result)) {
      cat(sprintf("No genes resolved for significant %s terms — plot skipped.\n", onto))
      next
    }

    out_file <- file.path(dir_violin, sprintf("GO_term_LFC_violin_ALL_%s_%s.png", onto, tag))
    ggsave(
      filename = out_file,
      plot     = result$plot,
      width    = 8,
      height   = max(6, 0.35 * result$n_terms + 2),
      dpi      = 300,
      limitsize = FALSE
    )
    cat(sprintf("Saved: %s (%d terms)\n", out_file, result$n_terms))
  }

  # 23i. Log2FC distribution violin + dots for BP "most-general parent" terms only
  bp_sig_all2 <- sig_results %>% dplyr::filter(Ontology == "BP")

  if (nrow(bp_sig_all2) > 0) {

    general_bp_ids2 <- filter_most_general_terms(bp_sig_all2$GO_ID, GOBPANCESTOR)

    bp_general_sig <- bp_sig_all2 %>%
      dplyr::filter(GO_ID %in% general_bp_ids2) %>%
      arrange(FDR)

    cat(sprintf("BP most-general terms for LFC plot: %d (out of %d significant BP terms)\n",
                nrow(bp_general_sig), nrow(bp_sig_all2)))

    result_bp_general <- make_all_terms_lfc_plot(
      onto_name       = "BP - most general terms only",
      sig_onto        = bp_general_sig,
      go2genes        = go2genes,
      gene_lfc_lookup = gene_lfc_lookup,
      contrast_label  = contrast_label,
      axis_label      = axis_label
    )

    if (!is.null(result_bp_general)) {
      out_file <- file.path(dir_violin, sprintf("GO_term_LFC_violin_BP_most_general_%s.png", tag))
      ggsave(
        filename  = out_file,
        plot      = result_bp_general$plot,
        width     = 8,
        height    = max(6, 0.35 * result_bp_general$n_terms + 2),
        dpi       = 300,
        limitsize = FALSE
      )
      cat(sprintf("Saved: %s (%d terms)\n", out_file, result_bp_general$n_terms))
    } else {
      cat("No genes resolved for most-general BP terms — plot skipped.\n")
    }

  } else {
    cat("No significant BP terms — most-general BP LFC violin plot skipped.\n")
  }

  cat(sprintf("\n Part B (camera GO enrichment) complete for %s.\n", contrast_label))

  invisible(list(res = res, camera_results = camera_results))
}

# 24. Run all contrasts
all_results <- list()
for (contrast_name in names(contrast_defs)) {
  all_results[[contrast_name]] <- run_contrast_analysis(contrast_name, contrast_defs[[contrast_name]])
}

# 25. Combined multi-contrast GO dotplot (one plot per ontology)
# x = contrast (LLvsDark, MLvsDark, HLvsDark, FLvsDark), y = GO term,
# colour = Direction, size = -log10(FDR)
cat("\n\n================ Building combined multi-contrast GO dotplots ================\n")

all_camera_results <- bind_rows(
  lapply(names(all_results), function(nm) {
    all_results[[nm]]$camera_results %>% mutate(Contrast = nm)
  })
)

combined_contrast_labels <- vapply(
  combined_plot_contrasts,
  function(nm) contrast_defs[[nm]]$label,
  character(1)
)

for (onto in ontologies) {

  result <- make_combined_dotplot(
    all_camera_results, onto,
    contrast_order  = combined_plot_contrasts,
    contrast_labels = combined_contrast_labels
  )

  if (is.null(result)) {
    cat(sprintf("No significant %s terms across %s — combined dotplot skipped.\n",
                onto, paste(combined_plot_contrasts, collapse = ", ")))
    next
  }

  out_file <- file.path(dir_combined, sprintf("camera_GO_dotplot_combined_%s_vsDark.png", onto))
  ggsave(
    filename  = out_file,
    plot      = result$plot,
    width     = 9,
    height    = max(5, 0.3 * result$n_terms + 2),
    dpi       = 300,
    limitsize = FALSE
  )
  cat(sprintf("Saved: %s (%d terms)\n", out_file, result$n_terms))
}

# 26. Combined multi-contrast BP dotplot, most-general terms only
# Same idea as 23f but across all contrasts at once (like the combined
# dotplots in 25): the most-general parent terms are determined from the
# union of significant BP terms across combined_plot_contrasts, then each
# surviving term is plotted across all of those contrasts.
cat("\n\n================ Building combined BP dotplot (most general terms only) ================\n")

bp_combined_sig_ids <- all_camera_results %>%
  dplyr::filter(Ontology == "BP", Contrast %in% combined_plot_contrasts, FDR < fdr_cutoff) %>%
  dplyr::pull(GO_ID) %>%
  unique()

if (length(bp_combined_sig_ids) == 0) {

  cat("No significant BP terms across combined contrasts — most-general combined BP dotplot skipped.\n")

} else {

  general_bp_combined_ids <- filter_most_general_terms(bp_combined_sig_ids, GOBPANCESTOR)

  cat(sprintf("Combined BP significant terms: %d total -> %d most-general (parents only, children dropped)\n",
              length(bp_combined_sig_ids), length(general_bp_combined_ids)))

  result_bp_general <- make_combined_dotplot(
    all_camera_results, "BP",
    contrast_order  = combined_plot_contrasts,
    contrast_labels = combined_contrast_labels,
    go_id_filter    = general_bp_combined_ids,
    title_label     = "BP - most general terms only"
  )

  if (is.null(result_bp_general)) {
    cat("No significant most-general BP terms across combined contrasts — dotplot skipped.\n")
  } else {
    out_file_bp_general <- file.path(dir_combined, "camera_GO_dotplot_combined_BP_most_general_vsDark.png")
    ggsave(
      filename  = out_file_bp_general,
      plot      = result_bp_general$plot,
      width     = 9,
      height    = max(5, 0.3 * result_bp_general$n_terms + 2),
      dpi       = 300,
      limitsize = FALSE
    )
    cat(sprintf("Saved: %s (%d terms)\n", out_file_bp_general, result_bp_general$n_terms))
  }
}

cat("\n Complete \n")
