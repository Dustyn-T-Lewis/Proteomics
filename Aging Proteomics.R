# ---- 0. Load Required Libraries ----
# Define and install required libraries
required_packages <- c(
  "readxl", "dplyr", "stringr", "tidyr", "httr", "jsonlite", "purrr",
  "tibble", "broom", "ReactomePA", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "ggplot2", "forcats", "fmsb", "BiocManager"
)

# Install missing libraries
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("ReactomePA", "clusterProfiler", "org.Hs.eg.db", "enrichplot")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

# Load libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# ---- 1. Load & Clean Data ----
# Load data from Excel, remove NA rows/columns, and extract gene/protein details
sarco <- read_excel("13_MASTER_SEER proteome analysis (4-4-23).xlsx", sheet = "MASTER_SARCO_data") %>%
  dplyr::select(where(~ any(!is.na(.)))) %>%
  filter(if_any(everything(), ~ !is.na(.))) %>%
  mutate(
    Gene_Symbol = str_remove(str_extract(Gene_protein, "^[^_;]+"), ";?GN=+$"),
    Accession   = str_extract(Gene_protein, "\\b[A-Z][0-9][A-Z0-9]{3}[0-9]\\b")
  )

# Calculate detection scores based on sample columns
sample_cols <- grep("^(EAA|PRE|POST|PPS)", names(sarco), value = TRUE)
sarco <- sarco %>% mutate(detection_score = rowSums(across(all_of(sample_cols)), na.rm = TRUE))

# ---- 2. Query UniProt for Missing Gene Symbols ----
# Define UniProt query function
query_uniprot <- function(accessions) {
  url <- "https://rest.uniprot.org/uniprotkb/search"
  q <- paste(accessions, collapse = " OR ")
  r <- GET(url, query = list(query = q, fields = "accession,gene_primary", format = "json", size = 500))
  if (status_code(r) == 200) {
    j <- httr::content(r, as = "parsed", type = "application/json")
    tibble(
      Accession = sapply(j$results, function(x) x$primaryAccession),
      Gene_Symbol_API = sapply(j$results, function(x)
        if (!is.null(x$genes[[1]]$geneName$value)) x$genes[[1]]$geneName$value else NA_character_)
    )
  } else stop("UniProt query failed")
}

# Identify accessions requiring a UniProt query
accessions_to_query <- sarco %>%
  group_by(Accession, Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(is.na(Gene_Symbol) | Gene_Symbol == "") %>%
  pull(Accession) %>% unique() %>% na.omit()

# Query UniProt and collect results
uniprot_results <- accessions_to_query %>%
  split(ceiling(seq_along(.) / 100)) %>%
  map_dfr(query_uniprot) %>%
  distinct(Accession, .keep_all = TRUE)

# ---- 3. Merge and Clean Gene Symbols ----
# Merge UniProt results and clean up gene symbols
sarco <- sarco %>%
  left_join(uniprot_results, by = "Accession") %>%
  mutate(
    Gene_Symbol = if_else(is.na(Gene_Symbol) | Gene_Symbol == "", Gene_Symbol_API, Gene_Symbol),
    Gene_Symbol = str_remove(Gene_Symbol, "^isoform_")
  ) %>%
  dplyr::select(-Gene_Symbol_API)

# ---- 4. Remove Duplicates by Detection Score ----
# Filter data to retain highest detection score per Gene/Accession
sarco_clean <- sarco %>%
  group_by(Accession, Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- 5. Convert to Long Format ----
# Reshape data for analysis
long_df <- sarco_clean %>%
  dplyr::select(-matches("^T-test_|^DELTA_|^CORREL_|^Young_|^detection_")) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(!str_detect(Sample, "_Avg|_StdDev")) %>%
  mutate(
    Timepoint = case_when(
      str_detect(Sample, "PRE") ~ "Pre",
      str_detect(Sample, "POST") ~ "Post",
      TRUE ~ NA_character_
    ),
    Group = case_when(
      str_detect(Sample, "EAA|PPS") ~ "Young",
      TRUE ~ "MA"
    )
  )

# Save cleaned data
write.csv(long_df, "Clean_Long2.csv", row.names = FALSE)
write.csv(sarco_clean, "Clean_Wide2.csv", row.names = FALSE)

# ---- 6. Summary Table: Mean Abundances ----
# Create summary table of mean abundances per group
mean_abundance_summary <- long_df %>%
  mutate(
    Grouping = case_when(
      Group == "Young" ~ "Young",
      Group == "MA" & Timepoint == "Pre"  ~ "Pre_MA",
      Group == "MA" & Timepoint == "Post" ~ "Post_MA",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Grouping)) %>%
  group_by(Gene_Symbol, Grouping) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Grouping, values_from = Mean_Abundance)

write.csv(mean_abundance_summary, "Group_Abundance.csv", row.names = FALSE)

# ---- 7. Compute log2FC and Per-Gene Stats ----
# Calculate log2 fold-change
fc_df <- mean_abundance_summary %>%
  mutate(
    log2FC_Pre_vs_Young   = log2(Pre_MA / Young),
    log2FC_Post_vs_Young  = log2(Post_MA / Young),
    log2FC_Post_vs_Pre_MA = log2(Post_MA / Pre_MA)
  )

# Perform statistical analysis with corrected logic
stat_df <- long_df %>%
  group_by(Gene_Symbol) %>%
  do({
    df <- .
    df1 <- df %>% filter((Group == "MA" & Timepoint == "Pre") | Group == "Young")
    p1  <- if (n_distinct(df1$Group) == 2) t.test(Abundance ~ Group, data = df1)$p.value else NA_real_
    df2 <- df %>% filter((Group == "MA" & Timepoint == "Post") | Group == "Young")
    p2  <- if (n_distinct(df2$Group) == 2) t.test(Abundance ~ Group, data = df2)$p.value else NA_real_
    df3 <- df %>% filter(Group == "MA", Timepoint %in% c("Pre", "Post"))
    p3  <- if (n_distinct(df3$Timepoint) == 2) t.test(Abundance ~ Timepoint, data = df3)$p.value else NA_real_
    tibble(
      p_Pre_vs_Young   = p1,
      p_Post_vs_Young  = p2,
      p_Post_vs_Pre_MA = p3
    )
  }) %>%
  ungroup() %>%
  mutate(
    adj_Pre_vs_Young   = p.adjust(p_Pre_vs_Young,   method = "fdr"),
    adj_Post_vs_Young  = p.adjust(p_Post_vs_Young,  method = "fdr"),
    adj_Post_vs_Pre_MA = p.adjust(p_Post_vs_Pre_MA, method = "fdr")
  )

gene_stats <- fc_df %>%
  left_join(stat_df, by = "Gene_Symbol")

# ---- 8. Pathway Enrichment Analysis (Aging in Muscle) ----
# Map gene symbols to Entrez IDs
entrez_map <- bitr(
  gene_stats$Gene_Symbol,
  fromType = "SYMBOL", toType = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# Prepare gene list for enrichment analysis (Aging in Muscle: Pre vs Young)
gene_list_Pre_vs_Young <- gene_stats %>%
  inner_join(entrez_map, by = c("Gene_Symbol" = "SYMBOL")) %>%
  arrange(desc(log2FC_Pre_vs_Young)) %>%
  { set_names(.$log2FC_Pre_vs_Young, .$ENTREZID) }

# Perform KEGG, Reactome, and GO enrichment analyses (Aging in Muscle: Pre vs Young)
kegg_Pre_vs_Young <- enrichKEGG(
  gene         = names(gene_list_Pre_vs_Young),
  organism     = "hsa",
  pvalueCutoff = 0.10
)

reactome_Pre_vs_Young <- enrichPathway(
  gene         = names(gene_list_Pre_vs_Young),
  organism     = "human",
  pvalueCutoff = 0.10
)

go_Pre_vs_Young <- enrichGO(
  gene         = names(gene_list_Pre_vs_Young),
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.10
)

# Store enrichment results in a named list
enrichment_results <- list(
  KEGG_Pre_vs_Young = kegg_Pre_vs_Young,
  Reactome_Pre_vs_Young = reactome_Pre_vs_Young,
  GO_Pre_vs_Young = go_Pre_vs_Young
)

# ---- 9. Visualize Enrichment Results (Aging in Muscle) ----
print(dotplot(kegg_Pre_vs_Young, showCategory = 15) +
        ggtitle("KEGG: Pre_MA vs Young (AGING)"))

print(dotplot(reactome_Pre_vs_Young, showCategory = 15) +
        ggtitle("Reactome: Pre_MA vs Young (AGING)"))

print(dotplot(go_Pre_vs_Young, showCategory = 15) +
        ggtitle("GO: Pre_MA vs Young (AGING)"))

# ---- 10. Facetted “Zoom‑In” Pathway Plots (p ≤ 0.10 only) ----

# Define pathway targets for aging and training
kegg_targets <- c(
  Glutathione = "Glutathione metabolism",
  Proteasome = "Proteasome",
  ER_Processing = "Protein processing in endoplasmic reticulum",
  Ubiquitin = "Ubiquitin mediated proteolysis"
)

reactome_targets <- c(
  DetoxROS = "Detoxification of Reactive Oxygen Species",
  ERAD = "ER-Associated Degradation (ERAD)",
  Folding = "Protein folding",
  Macroautophagy = "Macroautophagy"
)

go_targets <- c(
  OxidativeStress = "response to oxidative stress",
  ProteinFolding = "protein folding",
  Proteolysis = "proteasome-mediated ubiquitin-dependent protein catabolic process",
  UPR = "endoplasmic reticulum unfolded protein response"
)

# ---- 11.Function to create individual zoom-in plots with expanded thresholds
create_zoom_in_plot <- function(enrichment, pathway_name, target_description, contrast, db_name, max_genes = 10) {
  idx <- which(enrichment$Description == target_description)
  if (length(idx) == 0) idx <- grep(target_description, enrichment$Description, ignore.case = TRUE)
  if (length(idx) == 0) return(NULL)
  
  genes <- str_split(enrichment$geneID[idx[1]], "/")[[1]]
  symbols <- entrez_map$SYMBOL[match(genes, entrez_map$ENTREZID)]
  
  plot_data <- gene_stats %>%
    filter(Gene_Symbol %in% symbols, .data[[paste0("adj_", contrast)]] <= 0.20) %>%
    transmute(
      Gene_Symbol,
      value = .data[[paste0("log2FC_", contrast)]],
      star = case_when(
        .data[[paste0("adj_", contrast)]] < 0.001 ~ "***",
        .data[[paste0("adj_", contrast)]] < 0.01 ~ "**",
        .data[[paste0("adj_", contrast)]] < 0.05 ~ "*",
        .data[[paste0("adj_", contrast)]] < 0.10 ~ ".",
        TRUE ~ ""
      ),
      significance = case_when(
        .data[[paste0("adj_", contrast)]] < 0.05 ~ "Significant",
        .data[[paste0("adj_", contrast)]] < 0.10 ~ "Near Significant",
        .data[[paste0("adj_", contrast)]] < 0.20 ~ "Trending",
        TRUE ~ "Not Significant"
      )
    ) %>%
    arrange(desc(value)) %>%
    head(max_genes)
  
  if (nrow(plot_data) == 0) return(NULL)
  
  print(
    ggplot(plot_data, aes(x = reorder(Gene_Symbol, value), y = value, fill = significance)) +
      geom_col(show.legend = TRUE) +
      geom_text(aes(label = star), hjust = -0.3, size = 5, color = "black", fontface = "bold") +
      coord_flip() +
      scale_fill_manual(values = c(
        "Significant" = "#D55E00",
        "Near Significant" = "#F0E442",
        "Trending" = "#009E73",
        "Not Significant" = "#C0C0C0"
      )) +
      labs(
        title = paste(db_name, "Pathway:", pathway_name, "(", contrast, ")"),
        subtitle = paste(
          "Pathway:", target_description,
          "\nSignificance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10 (Near Significant), Marginal p < 0.20"
        ),
        x = "Gene Symbol",
        y = paste0("log₂ Fold Change (", contrast, ")"),
        caption = "Bars represent log2 fold change for each gene in the pathway.\nColored by significance threshold."
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 14, margin = margin(b = 10)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.caption = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
      )
  )
}

# ---- 12. Generate plots for each database and contrast
contrasts <- c("Pre_vs_Young", "Post_vs_Pre_MA") # Aging and training contrasts

# KEGG
for (contrast in contrasts) {
  for (pathway_name in names(kegg_targets)) {
    enrichment <- enrichment_results[[paste0("KEGG_", contrast)]]
    create_zoom_in_plot(enrichment, pathway_name, kegg_targets[[pathway_name]], contrast, "KEGG", max_genes = 20)
  }
}

# Reactome
for (contrast in contrasts) {
  for (pathway_name in names(reactome_targets)) {
    enrichment <- enrichment_results[[paste0("Reactome_", contrast)]]
    create_zoom_in_plot(enrichment, pathway_name, reactome_targets[[pathway_name]], contrast, "Reactome", max_genes = 20)
  }
}

# GO
for (contrast in contrasts) {
  for (pathway_name in names(go_targets)) {
    enrichment <- enrichment_results[[paste0("GO_", contrast)]]
    create_zoom_in_plot(enrichment, pathway_name, go_targets[[pathway_name]], contrast, "GO", max_genes = 20)
  }
}

