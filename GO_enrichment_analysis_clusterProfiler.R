#!/usr/bin/env Rscript
# ============================================================================
# GO Term Enrichment Analysis — Chlamydomonas reinhardtii, FL vs HL clusters
# Using clusterProfiler::enricher() with GO hierarchy propagation
# Background (universe) = whole genome (all genes in the dataset)
# ============================================================================

# ---- A. Packages -----------------------------------------------------------
suppressMessages({
  library(clusterProfiler)
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
out_dir    <- "GO_enrichment_results"
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
# just genes with GO annotation, as requested.

# whole_genome <- unique(df$gene_name)
# cat("Whole genome size (universe):", length(whole_genome), "genes\n")

# ---- B2. Build direct gene -> GO term mapping (only annotated genes) -------
gene2go_raw <- df %>%
  filter(!is.na(go_terms), str_trim(go_terms) != "") %>%
  select(gene_name, go_terms)

# Background - annotated genes
whole_genome <- gene2go_raw$gene_name %>% unique()

# GO Term column is space-separated -> split into individual terms
gene2go_direct <- setNames(
  str_split(gene2go_raw$go_terms, "\\s+"),
  gene2go_raw$gene_name
)
gene2go_direct <- lapply(gene2go_direct, function(x) unique(x[x != ""]))

cat("Genes with direct GO annotation:", length(gene2go_direct), "\n")

# ---- B3. Map each direct GO term to its ontology (BP/MF/CC) -----------------
all_direct_terms <- unique(unlist(gene2go_direct, use.names = FALSE))
cat("Unique direct GO terms:", length(all_direct_terms), "\n")

# Only keep terms that actually exist in GO.db (drops obsolete/typo'd IDs)
valid_terms <- all_direct_terms[all_direct_terms %in% keys(GO.db, keytype = "GOID")]
dropped <- setdiff(all_direct_terms, valid_terms)
if (length(dropped) > 0) {
  cat("Dropping", length(dropped), "GO terms not found in GO.db (obsolete/invalid)\n")
}

go_info <- AnnotationDbi::select(GO.db, keys = valid_terms,
                                  keytype = "GOID",
                                  columns = c("ONTOLOGY", "TERM"))
id2onto <- setNames(go_info$ONTOLOGY, go_info$GOID)
id2name <- setNames(go_info$TERM, go_info$GOID)

# Restrict gene2go_direct to valid terms only
gene2go_direct <- lapply(gene2go_direct, function(x) intersect(x, valid_terms))
gene2go_direct <- gene2go_direct[lengths(gene2go_direct) > 0]

# ============================================================================
# B4. Propagate GO hierarchy: gene -> direct + all ancestor terms
# ============================================================================
# Expand ancestors ONCE per unique GO term (not per gene), then do fast
# list lookups for each gene. Avoids thousands of redundant mget() calls.
cat("Propagating GO hierarchy (BP/MF/CC ancestors)...\n")
ancestor_dbs <- list(BP = GOBPANCESTOR, MF = GOMFANCESTOR, CC = GOCCANCESTOR)

# Build term -> (self + all ancestors) lookup, one ontology at a time
term2full <- list()
for (onto in names(ancestor_dbs)) {
  terms_onto <- names(id2onto)[id2onto == onto]
  if (length(terms_onto) == 0) next
  db       <- ancestor_dbs[[onto]]
  mapped   <- terms_onto[terms_onto %in% mappedkeys(db)]
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
cat("Hierarchy propagation done.\n")

# ---- B5. Recompute ontology/name for the FULL propagated term set ----------
# (ancestor terms weren't necessarily in the original valid_terms/go_info set)
all_full_terms <- unique(unlist(gene2go_full, use.names = FALSE))
new_terms <- setdiff(all_full_terms, names(id2onto))
if (length(new_terms) > 0) {
  go_info_new <- AnnotationDbi::select(GO.db, keys = new_terms,
                                        keytype = "GOID",
                                        columns = c("ONTOLOGY", "TERM"))
  id2onto <- c(id2onto, setNames(go_info_new$ONTOLOGY, go_info_new$GOID))
  id2name <- c(id2name, setNames(go_info_new$TERM,     go_info_new$GOID))
}

# ---- B6. Build TERM2GENE / TERM2NAME tables per ontology -------------------
# Long-format gene -> GO term table (post-propagation)
gene2go_long <- tibble(
  gene = rep(names(gene2go_full), lengths(gene2go_full)),
  term = unlist(gene2go_full, use.names = FALSE)
) %>%
  mutate(ontology = id2onto[term]) %>%
  filter(!is.na(ontology))  # safety: drop any term we couldn't classify

build_term2gene <- function(onto) {
  gene2go_long %>% filter(ontology == onto) %>% select(term, gene)
}
build_term2name <- function(onto) {
  terms_onto <- unique(gene2go_long$term[gene2go_long$ontology == onto])
  tibble(term = terms_onto, name = id2name[terms_onto])
}

term2gene_list <- setNames(lapply(c("BP", "MF", "CC"), build_term2gene), c("BP", "MF", "CC"))
term2name_list <- setNames(lapply(c("BP", "MF", "CC"), build_term2name), c("BP", "MF", "CC"))

for (onto in c("BP", "MF", "CC")) {
  cat(sprintf("  %s: %d gene-term associations, %d unique terms\n",
              onto, nrow(term2gene_list[[onto]]), nrow(term2name_list[[onto]])))
}

# ============================================================================
# C. Run enrichment per cluster x ontology
# ============================================================================
clusters <- df %>% filter(!is.na(cluster), str_trim(cluster) != "") %>%
  pull(cluster) %>% unique() %>% sort()
cat("Clusters found:", paste(clusters, collapse = ", "), "\n")

run_enrichment <- function(gene_list, universe, term2gene, term2name, onto, cluster_name) {
  # enricher requires at least a couple genes overlapping the term2gene universe
  tryCatch({
    ego <- enricher(
      gene          = gene_list,
      universe      = universe,
      TERM2GENE     = term2gene,
      TERM2NAME     = term2name,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,     # keep all results, filter later; adjust as desired
      qvalueCutoff  = 1,
      minGSSize     = 5,
      maxGSSize     = 500
    )
    if (is.null(ego) || nrow(ego@result) == 0) {
      cat(sprintf("    [%s / %s] no enrichment result\n", cluster_name, onto))
      return(NULL)
    }
    res <- ego@result %>%
      mutate(Cluster = cluster_name, Ontology = onto) %>%
      relocate(Cluster, Ontology)
    res
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
    res <- run_enrichment(
      gene_list = gene_list,
      universe  = whole_genome,           # whole genome as background
      term2gene = term2gene_list[[onto]],
      term2name = term2name_list[[onto]],
      onto      = onto,
      cluster_name = cl
    )
    if (!is.null(res)) {
      all_results[[paste(cl, onto, sep = "_")]] <- res
    }
  }
}

# ---- D. Combine, apply significance threshold, and export -----------------
combined <- bind_rows(all_results)

# Multiple-testing correction strategy: BH correction applied per (Cluster x
# Ontology) test above is already the p.adjust column from enricher(). If you
# prefer BH correction pooled across all clusters/ontologies simultaneously
# (often more defensible for a single manuscript-wide claim), recompute below:
combined <- combined %>%
  mutate(p.adjust.pooled = p.adjust(pvalue, method = "BH"))

sig_results <- combined %>% filter(p.adjust.pooled < 0.05)

write_csv(combined,     file.path(out_dir, "GO_enrichment_all_results.csv"))
write_csv(sig_results,  file.path(out_dir, "GO_enrichment_significant_results.csv"))

cat("\n============================================================\n")
cat("Done.\n")
cat("All results:        ", file.path(out_dir, "GO_enrichment_all_results.csv"), "\n")
cat("Significant (p.adj < 0.05):", file.path(out_dir, "GO_enrichment_significant_results.csv"), "\n")
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
    filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    slice_head(n = top_n) %>%
    mutate(
      # GeneRatio comes back from enricher() as "k/n" string -> convert to numeric
      GeneRatio = sapply(GeneRatio, function(x) eval(parse(text = x))),
      # Truncate long GO term names so labels stay legible on the plot
      Description = if_else(nchar(Description) > 50, paste0(str_sub(Description, 1, 47), "..."), Description),
      Description = fct_reorder(Description, GeneRatio)
    )

  if (nrow(d) == 0) {
    cat(sprintf("  [%s / %s] no terms to plot, skipping\n", cluster_name, onto))
    return(NULL)
  }

  p <- ggplot(d, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust,pooled)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    scale_size_continuous(name = "Gene\ncount") +
    labs(
      title = sprintf("GO Enrichment (%s) - Cluster %s", onto, cluster_name),
      x = "Gene Ratio",
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
