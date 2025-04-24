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

# Function to create individual zoom-in plots with expanded thresholds
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
        .data[[paste0("adj_", contrast)]] < 0.20 ~ "Marginally Significant",
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
        "Marginally Significant" = "#009E73",
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

# Generate plots for each database and contrast
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

# ---- 11. Individual Pathway Radar and Bar Plots with Expanded Categories ----

# Define expanded pathway categories
categories <- list(
  Redox = list(
    Redox_Stress = c("NOX1", "NOX2", "NOX4", "CYBA", "CYBB", "DUOX1", "DUOX2", "XDH", "AIFM1", "AIFM2"),
    Antioxidant_Defenses = c("SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX2", "GPX3", "GPX4",
                             "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
                             "TXN", "TXNRD1", "TXNRD2", "GSR", "GCLC", "GCLM", "HMOX1"),
    Oxidative_Damage_Signalling = c("NFE2L2", "KEAP1", "NQO1", "HMOX1", "TXNIP", "PARP1", "ATM",
                                    "ATR", "TP53", "NFKB1", "RELA")
  ),
  Synthesis = list(
    Protein_Synthesis_Signalling = c("MTOR", "RPTOR", "RICTOR", "AKT1", "AKT2", "TSC1", "TSC2",
                                     "RHEB", "PRKAA1", "PRKAA2", "RPS6KB1", "RPS6KB2",
                                     "EIF4EBP1", "EIF4E"),
    Ribosome_Biogenesis = c(paste0("RPL", c("P0", "P1", "P2", "3", "4", "7", "11")),
                            paste0("RPS", c("3", "6", "8", "19")),
                            "MRPL12", "NOP14", "NOP56", "UTP15", "FBL", "PES1", "BOP1", "WDR12"),
    ER_Stress = c("HSPA5", "HSP90B1", "CANX", "CALR", "PDIA1", "PDIA3", "PDIA4",
                  "PDIA6", "DDIT3", "GADD34", "ATF6", "EIF2AK3", "ERN1", "ERN2",
                  "SEC61A1", "DERL1"),
    UPR_and_Chaperones = c("XBP1", "XBP1s", "ATF4", "DDIT3", "HSPA5", "HSP90AA1", "HSP90AB1",
                           "DNAJA1", "DNAJB1", "DNAJB9", "DNAJC3", "BAG3", "GRP94")
  ),
  Proteolysis = list(
    Calpains = c("CAPN1", "CAPN2", "CAPN3", "CAPN5", "CAPN10", "CAPNS1", "CAST"),
    Caspases = paste0("CASP", c(1, 3, 4, 5, 6, 7, 8, 9, 10, 12)),
    Autophagy = c("ULK1", "BECN1", "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1",
                  "MAP1LC3A", "MAP1LC3B", "SQSTM1"),
    Ubiquitin_Proteasome = c(paste0("PSMA", 1:7), paste0("PSMB", 1:8),
                             paste0("PSMC", 1:6), paste0("PSMD", 1:14),
                             "UBB", "UBC", "UBA1", "UBE2D2", "UBE2L3", "UBE3A")
  ),
  Mito = list(
    Fusion = c("MFN1", "MFN2", "OPA1", "OPA3"),
    Fission = c("DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2"),
    Mitophagy = c("PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "BCL2L13"),
    Biogenesis = c("PPARGC1A", "PPARGC1B", "NRF1", "NRF2", "TFAM", "TFB2M")
  )
)

# Define groups and colors
groups <- c("Young", "Pre", "Post")
colors <- c("#1b9e77", "#d95f02", "#7570b3")
group_order <- c("Young", "Pre", "Post")

# Helper function to calculate data for each pathway
calculate_pathway_data <- function(category, data) {
  abundance_matrix <- sapply(groups, function(group) {
    subset_data <- data %>%
      filter(
        (group == "Young" & Group == "Young") |
          (group == "Pre" & Group == "MA" & Timepoint == "Pre") |
          (group == "Post" & Group == "MA" & Timepoint == "Post")
      )
    sapply(category, function(genes) mean(subset_data$Abundance[subset_data$Gene_Symbol %in% genes], na.rm = TRUE))
  })
  list(
    normalized = abundance_matrix / abundance_matrix[, "Young"],
    raw = abundance_matrix
  )
}

# Prepare data for pathway categories
data_redox <- calculate_pathway_data(categories$Redox, long_df)
data_synthesis <- calculate_pathway_data(categories$Synthesis, long_df)
data_proteolysis <- calculate_pathway_data(categories$Proteolysis, long_df)
data_mito <- calculate_pathway_data(categories$Mito, long_df)

# Generate radar plots
plot_radar <- function(data, pathway_name) {
  radar_data <- as.data.frame(t(data$normalized)) %>%
    rownames_to_column("Group") %>%
    pivot_longer(-Group, names_to = "Module", values_to = "Abundance")
  
  radar_plot <- ggplot(radar_data, aes(x = Module, y = Abundance, color = Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    coord_polar() +
    labs(
      title = paste(pathway_name, "Radar Plot"),
      x = "Pathway Module",
      y = "Normalized Abundance",
      color = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  print(radar_plot)
}

# Generate bar plots
plot_bar <- function(data, pathway_name) {
  bar_data <- as.data.frame(t(data$raw)) %>%
    rownames_to_column("Group") %>%
    pivot_longer(-Group, names_to = "Module", values_to = "Abundance") %>%
    mutate(Group = factor(Group, levels = group_order))
  
  bar_plot <- ggplot(bar_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_col(position = "dodge", color = "black") +
    facet_wrap(~ Module, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste(pathway_name, "Pathway Abundance Bar Plot"),
      x = "Group",
      y = "Mean Abundance",
      fill = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  print(bar_plot)
}

# Create and print radar and bar plots for each pathway
plot_radar(data_redox, "Redox")
plot_bar(data_redox, "Redox")
plot_radar(data_synthesis, "Synthesis")
plot_bar(data_synthesis, "Synthesis")
plot_radar(data_proteolysis, "Proteolysis")
plot_bar(data_proteolysis, "Proteolysis")
plot_radar(data_mito, "Mito")
plot_bar(data_mito, "Mito")

# ---- 12. Save Results ----
# Save enrichment analysis results
write.csv(kegg_Pre_vs_Young@result, "KEGG_Pre_vs_Young_Results.csv", row.names = FALSE)
write.csv(reactome_Pre_vs_Young@result, "Reactome_Pre_vs_Young_Results.csv", row.names = FALSE)
write.csv(go_Pre_vs_Young@result, "GO_Pre_vs_Young_Results.csv", row.names = FALSE)

# Save processed pathway data
saveRDS(data_redox, "Pathway_Redox_Data.rds")
saveRDS(data_synthesis, "Pathway_Synthesis_Data.rds")
saveRDS(data_proteolysis, "Pathway_Proteolysis_Data.rds")
saveRDS(data_mito, "Pathway_Mito_Data.rds")

# Save radar and bar plots as PDFs
pdf("Redox_Pathway_Plots.pdf")
plot_radar(data_redox, "Redox")
plot_bar(data_redox, "Redox")
dev.off()

pdf("Synthesis_Pathway_Plots.pdf")
plot_radar(data_synthesis, "Synthesis")
plot_bar(data_synthesis, "Synthesis")
dev.off()

pdf("Proteolysis_Pathway_Plots.pdf")
plot_radar(data_proteolysis, "Proteolysis")
plot_bar(data_proteolysis, "Proteolysis")
dev.off()

pdf("Mito_Pathway_Plots.pdf")
plot_radar(data_mito, "Mito")
plot_bar(data_mito, "Mito")
dev.off()