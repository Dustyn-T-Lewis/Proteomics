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

# ---- 6. Define Pathway Categories ----
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