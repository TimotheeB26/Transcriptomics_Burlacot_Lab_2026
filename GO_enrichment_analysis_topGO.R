#!/usr/bin/env Rscript
# ============================================================================
# GO Term Enrichment Analysis — Chlamydomonas reinhardtii, FL vs HL clusters
# Using topGO (weight01 algorithm) — GO hierarchy propagation handled natively
# by topGO's algorithm via the GO graph structure.
# Background (universe) = whole genome (all genes in the dataset)
# ============================================================================

# ---- A. Packages -----------------------------------------------------------
suppressMessages({
  library(topGO)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# ---- B. Load data -----------------------------------------------------------
input_file <- "Dataset_GO_Enr.csv"   # adjust path as needed
out_dir    <- "GO_enrichment_results_topGO"
dir.create(out_dir, showWarnings = FALSE)

df <- read_csv(input_file, show_col_types = FALSE)

# Expect columns: `Gene Name`, `Gene ID`, `GO Term`, `Predalgo`, `Cluster`
# NOTE: per your dataset, `Gene ID` holds the Phytozome/JGI v5 identifiers
# (CHLRE_*v5) and `Gene Name` holds the Cre-locus identifiers (Cre##.g######).
# We use `Gene Name` (Cre locus) as the gene identifier throughout since it's
# the more standard/readable ID, but you can switch to `Gene ID` if preferred.

df <- df %>%
  rename(gene_name = `Gene Name`, gene_id = `Gene ID`,
         go_terms  = `GO Term`, cluster = Cluster)

# ---- B1. Whole-genome universe ---------------------------------------------
# The universe/background = every gene in the dataset (whole genome), not
# just genes with GO annotation, as requested. topGO's allGenes vector must
# cover this full universe; genes with no GO annotation simply contribute no
# information to the tests but stay in the background.

# whole_genome <- unique(df$gene_name)
# cat("Whole genome size (universe):", length(whole_genome), "genes\n")

# ---- B2. Build direct gene -> GO term mapping (only annotated genes) -------
# topGO's weight01 algorithm propagates the GO hierarchy internally via the
# GO DAG (using the annotation you give it as the "direct" annotation, then
# walking up parent/child relationships itself during testing) — so we only
# need to supply DIRECT annotations here, not pre-propagated ones.
gene2go_raw <- df %>%
  filter(!is.na(go_terms), str_trim(go_terms) != "") %>%
  select(gene_name, go_terms)

# Background - annotated genes
whole_genome <- gene2go_raw$gene_name %>% unique()

gene2go_direct <- setNames(
  str_split(gene2go_raw$go_terms, "\\s+"),
  gene2go_raw$gene_name
)
gene2go_direct <- lapply(gene2go_direct, function(x) unique(x[x != ""]))

# Keep only GO IDs that exist in GO.db (drops obsolete/typo'd IDs)
all_direct_terms <- unique(unlist(gene2go_direct, use.names = FALSE))
valid_terms <- all_direct_terms[all_direct_terms %in% keys(GO.db, keytype = "GOID")]
dropped <- setdiff(all_direct_terms, valid_terms)
if (length(dropped) > 0) {
  cat("Dropping", length(dropped), "GO terms not found in GO.db (obsolete/invalid)\n")
}
gene2go_direct <- lapply(gene2go_direct, function(x) intersect(x, valid_terms))
gene2go_direct <- gene2go_direct[lengths(gene2go_direct) > 0]

cat("Genes with direct GO annotation:", length(gene2go_direct), "\n")

# ---- B3. Assign ontology (BP/MF/CC) to each direct term --------------------
go_info <- AnnotationDbi::select(GO.db, keys = valid_terms,
                                  keytype = "GOID",
                                  columns = c("ONTOLOGY", "TERM"))
id2onto <- setNames(go_info$ONTOLOGY, go_info$GOID)

# Build per-ontology gene2GO lists (topGO needs one list per ontology, since
# a gene's BP/MF/CC annotations are tested against separate GO graphs)
build_gene2go_onto <- function(onto) {
  terms_onto <- names(id2onto)[id2onto == onto]
  g2g <- lapply(gene2go_direct, function(x) intersect(x, terms_onto))
  g2g[lengths(g2g) > 0]
}
gene2go_by_onto <- setNames(lapply(c("BP", "MF", "CC"), build_gene2go_onto),
                             c("BP", "MF", "CC"))

for (onto in c("BP", "MF", "CC")) {
  cat(sprintf("  %s: %d genes with direct annotation\n", onto, length(gene2go_by_onto[[onto]])))
}

# ============================================================================
# C. Run topGO (weight01, Fisher's exact test) per cluster x ontology
# ============================================================================
clusters <- df %>% filter(!is.na(cluster), str_trim(cluster) != "") %>%
  pull(cluster) %>% unique() %>% sort()
cat("Clusters found:", paste(clusters, collapse = ", "), "\n")

run_topgo <- function(gene_list, universe, gene2go_onto, onto, cluster_name) {
  # allGenes: named factor over the WHOLE GENOME (universe), 1 = "interesting"
  # (in this cluster), 0 = not. This is how topGO encodes the background.
  all_genes <- factor(as.integer(universe %in% gene_list))
  cluster_total_annotated <- sum(names(gene2go_onto) %in% gene_list)
  names(all_genes) <- universe

  tryCatch({
    GOdata <- new(
      "topGOdata",
      ontology      = onto,
      allGenes      = all_genes,
      geneSelectionFun = function(x) x == 1,
      annot         = annFUN.gene2GO,
      gene2GO       = gene2go_onto,
      nodeSize      = 5      # ignore GO terms with fewer than 5 annotated genes
    )

    result_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

    n_terms <- length(usedGO(GOdata))
    res_table <- GenTable(
      GOdata,
      weight01 = result_weight01,
      orderBy  = "weight01",
      topNodes = n_terms
    )

    res_table <- res_table %>%
      as_tibble() %>%
      mutate(
        weight01 = as.numeric(gsub("< ", "", weight01)),
        Cluster  = cluster_name,
        Ontology = onto,
        ClusterTotalAnnotated = cluster_total_annotated   # <-- add this line
      ) %>%
      relocate(Cluster, Ontology)

    res_table
  }, error = function(e) {
    cat(sprintf("    [%s / %s] error: %s\n", cluster_name, onto, conditionMessage(e)))
    NULL
  })
}

all_results <- list()

for (cl in clusters) {
  gene_list <- df %>% filter(cluster == cl) %>% pull(gene_name) %>% unique()
  cat(sprintf("\nCluster %s: %d genes\n", cl, length(gene_list)))

  for (onto in c("BP", "MF", "CC")) {
    res <- run_topgo(
      gene_list    = gene_list,
      universe     = whole_genome,
      gene2go_onto = gene2go_by_onto[[onto]],
      onto         = onto,
      cluster_name = cl
    )
    if (!is.null(res)) {
      all_results[[paste(cl, onto, sep = "_")]] <- res
    }
  }
}

# ---- D. Combine, apply BH correction, and export ---------------------------
combined <- bind_rows(all_results)

# topGO's weight01 p-values are NOT independent test p-values in the usual
# sense (the algorithm already adjusts significance by removing genes
# assigned to significant child terms before testing parents), so applying
# a second BH correction on top is debated in the literature. We provide both
# the raw weight01 p-value and a BH-adjusted version (pooled across clusters
# WITHIN each ontology, i.e. one BH correction per BP/MF/CC) so you can choose
# your preferred convention, and note this caveat prominently in any
# report/manuscript methods section.
combined <- combined %>%
  group_by(Ontology) %>%
  mutate(weight01.padj.pooled = p.adjust(weight01, method = "BH")) %>%
  ungroup()

write_csv(combined, file.path(out_dir, "GO_enrichment_topGO_all_results.csv"))

sig_results <- combined %>% filter(weight01 < 0.05)
write_csv(sig_results, file.path(out_dir, "GO_enrichment_topGO_significant_raw.csv"))

sig_results_padj <- combined %>% filter(weight01.padj.pooled < 0.05)
write_csv(sig_results_padj, file.path(out_dir, "GO_enrichment_topGO_significant_padj.csv"))

cat("\n============================================================\n")
cat("Done.\n")
cat("All results:                 ", file.path(out_dir, "GO_enrichment_topGO_all_results.csv"), "\n")
cat("Significant (raw p < 0.05):  ", file.path(out_dir, "GO_enrichment_topGO_significant_raw.csv"), "\n")
cat("Significant (BH p.adj<0.05): ", file.path(out_dir, "GO_enrichment_topGO_significant_padj.csv"), "\n")
cat("============================================================\n")

# ============================================================================
# E. Dotplots — one per cluster x ontology
# ============================================================================
suppressMessages({
  library(ggplot2)
  library(forcats)
})

plot_dir <- file.path(out_dir, "dotplots")
dir.create(plot_dir, showWarnings = FALSE)

TOP_N <- 15  # number of top GO terms to show per dotplot

make_dotplot <- function(cluster_name, onto, data, top_n = TOP_N) {
  d <- data %>%
    filter(Cluster == cluster_name, Ontology == onto) %>%
    filter(!is.na(weight01)) %>%
    arrange(weight01) %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatio = Significant / ClusterTotalAnnotated,
      # Truncate long GO term names so labels stay legible on the plot
      Term = if_else(nchar(Term) > 50, paste0(str_sub(Term, 1, 47), "..."), Term),
      Term = fct_reorder(Term, GeneRatio)
    )

  if (nrow(d) == 0) {
    cat(sprintf("  [%s / %s] no terms to plot, skipping\n", cluster_name, onto))
    return(NULL)
  }

  p <- ggplot(d, aes(x = GeneRatio, y = Term, size = Significant, color = weight01)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "weight01\np-value") +
    scale_size_continuous(name = "Gene\ncount") +
    labs(
      title = sprintf("GO Enrichment (%s) - Cluster %s", onto, cluster_name),
      x = "Gene Ratio (Significant / Cluster size)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )

  fname <- file.path(plot_dir, sprintf("dotplot_%s_%s.png", cluster_name, onto))
  ggsave(fname, plot = p, width = 8, height = max(4, 0.35 * nrow(d) + 1.5), dpi = 300)
  cat(sprintf("  [%s / %s] saved: %s\n", cluster_name, onto, fname))
  p
}

cat("\nGenerating dotplots...\n")
for (cl in clusters) {
  for (onto in c("BP", "MF", "CC")) {
    make_dotplot(cl, onto, combined)
  }
}
cat("Dotplots saved in:", plot_dir, "\n")
