#### Proteomic Data Processing Pipeline (Concise & Complete) ####

# ---- 0. Load Required Libraries ----
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(httr)
library(jsonlite)
library(purrr)
library(tibble)
library(broom)
if(!requireNamespace("ReactomePA", quietly=TRUE)) BiocManager::install("ReactomePA")
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(forcats)
library(fmsb)

# ---- 1. Load & Clean Data ----
sarco <- read_excel("13_MASTER_SEER proteome analysis (4-4-23).xlsx", sheet = "MASTER_SARCO_data") %>%
  dplyr::select(where(~ any(!is.na(.)))) %>%  # Remove all-NA columns
  filter(if_any(everything(), ~ !is.na(.))) %>%  # Remove all-NA rows
  mutate(
    Gene_Symbol = str_remove(str_extract(Gene_protein, "^[^_;]+"), ";?GN=+$"),
    Accession   = str_extract(Gene_protein, "\\b[A-Z][0-9][A-Z0-9]{3}[0-9]\\b")
  )

sample_cols <- grep("^(EAA|PRE|POST|PPS)", names(sarco), value = TRUE)
sarco <- sarco %>% mutate(detection_score = rowSums(across(all_of(sample_cols)), na.rm = TRUE))

# ---- 2. Query UniProt for Missing Gene Symbols ----
accessions_to_query <- sarco %>%
  group_by(Accession, Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(is.na(Gene_Symbol) | Gene_Symbol == "") %>%
  pull(Accession) %>% unique() %>% na.omit()

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

uniprot_results <- accessions_to_query %>%
  split(ceiling(seq_along(.) / 100)) %>%
  map_dfr(query_uniprot) %>%
  distinct(Accession, .keep_all = TRUE)

# ---- 3. Merge and Clean Gene Symbols ----
sarco <- sarco %>%
  left_join(uniprot_results, by = "Accession") %>%
  mutate(
    Gene_Symbol = if_else(is.na(Gene_Symbol) | Gene_Symbol == "", Gene_Symbol_API, Gene_Symbol),
    Gene_Symbol = str_remove(Gene_Symbol, "^isoform_")
  ) %>%
  dplyr::select(-Gene_Symbol_API)

# ---- 4. Remove Duplicates by Detection Score ----
sarco_clean <- sarco %>%
  group_by(Accession, Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- 5. Convert to Long Format ----
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

write.csv(long_df, "Clean_Long2.csv", row.names = FALSE)
write.csv(sarco_clean, "Clean_Wide2.csv", row.names = FALSE)

# ---- 6. Summary Table: Mean Abundances ----
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
fc_df <- mean_abundance_summary %>%
  mutate(
    log2FC_Pre_vs_Young    = log2(Pre_MA  / Young),
    log2FC_Post_vs_Young   = log2(Post_MA / Young),
    log2FC_Post_vs_Pre_MA  = log2(Post_MA / Pre_MA)
  )

stat_df <- long_df %>%
  group_by(Gene_Symbol) %>%
  do({
    df <- .
    df1 <- df %>% filter((Group=="MA" & Timepoint=="Pre") | Group=="Young")
    p1  <- if (n_distinct(df1$Group)==2) t.test(Abundance ~ Group, data = df1)$p.value else NA_real_
    df2 <- df %>% filter((Group=="MA" & Timepoint=="Post") | Group=="Young")
    p2  <- if (n_distinct(df2$Group)==2) t.test(Abundance ~ Group, data = df2)$p.value else NA_real_
    df3 <- df %>% filter(Group=="MA", Timepoint %in% c("Pre","Post"))
    p3  <- if (n_distinct(df3$Timepoint)==2) t.test(Abundance ~ Timepoint, data = df3)$p.value else NA_real_
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

# ---- 8. Pathway Enrichment Analysis ----
entrez_map <- bitr(
  gene_stats$Gene_Symbol,
  fromType = "SYMBOL", toType = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_list_Pre_vs_Young <- gene_stats %>%
  inner_join(entrez_map, by = c("Gene_Symbol" = "SYMBOL")) %>%
  arrange(desc(log2FC_Pre_vs_Young)) %>%
  { set_names(.$log2FC_Pre_vs_Young, .$ENTREZID) }

gene_list_Post_vs_Pre_MA <- gene_stats %>%
  inner_join(entrez_map, by = c("Gene_Symbol" = "SYMBOL")) %>%
  arrange(desc(log2FC_Post_vs_Pre_MA)) %>%
  { set_names(.$log2FC_Post_vs_Pre_MA, .$ENTREZID) }

kegg_Pre_vs_Young <- enrichKEGG(
  gene         = names(gene_list_Pre_vs_Young),
  organism     = "hsa",
  pvalueCutoff = 0.10
)
kegg_Post_vs_Pre_MA <- enrichKEGG(
  gene         = names(gene_list_Post_vs_Pre_MA),
  organism     = "hsa",
  pvalueCutoff = 0.10
)
reactome_Pre_vs_Young <- enrichPathway(
  gene         = names(gene_list_Pre_vs_Young),
  organism     = "human",
  pvalueCutoff = 0.10
)
reactome_Post_vs_Pre_MA <- enrichPathway(
  gene         = names(gene_list_Post_vs_Pre_MA),
  organism     = "human",
  pvalueCutoff = 0.10
)
go_Pre_vs_Young <- enrichGO(
  gene         = names(gene_list_Pre_vs_Young),
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.10
)
go_Post_vs_Pre_MA <- enrichGO(
  gene         = names(gene_list_Post_vs_Pre_MA),
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.10
)

# ---- 9. Pathway Dotplots ----
dotplot(kegg_Pre_vs_Young,     showCategory=15) + ggtitle("KEGG: Pre_MA vs Young (AGING)")
dotplot(reactome_Pre_vs_Young, showCategory=15) + ggtitle("Reactome: Pre_MA vs Young (AGING)")
dotplot(go_Pre_vs_Young,       showCategory=15) + ggtitle("GO: Pre_MA vs Young (AGING)")
dotplot(kegg_Post_vs_Pre_MA,     showCategory=15) + ggtitle("KEGG: MA Post vs MA Pre (TRAINING)")
dotplot(reactome_Post_vs_Pre_MA, showCategory=15) + ggtitle("Reactome: MA Post vs MA Pre (TRAINING)")
dotplot(go_Post_vs_Pre_MA,       showCategory=15) + ggtitle("GO: MA Post vs MA Pre (TRAINING)")

# 10. Facetted “zoom‑in” plots by database & contrast (p ≤ 0.10 only) 
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

# → Define the pathways to zoom into
kegg_targets <- c(
  Glutathione   = "Glutathione metabolism",
  Proteasome    = "Proteasome",
  ER_Processing = "Protein processing in endoplasmic reticulum",
  Ubiquitin     = "Ubiquitin mediated proteolysis"
)
reactome_targets <- c(
  DetoxROS       = "Detoxification of Reactive Oxygen Species",
  ERAD           = "ER-Associated Degradation (ERAD)",
  Folding        = "Protein folding",
  Macroautophagy = "Macroautophagy"
)
go_targets <- c(
  OxidativeStress = "response to oxidative stress",
  ProteinFolding  = "protein folding",
  Proteolysis     = "proteasome-mediated ubiquitin-dependent protein catabolic process",
  UPR             = "endoplasmic reticulum unfolded protein response"
)

# helper to pull SYMBOLs out of any enrichResult
get_syms <- function(enr, desc) {
  idx <- which(enr$Description == desc)
  if (length(idx)==0) idx <- grep(desc, enr$Description, ignore.case=TRUE)
  if (length(idx)==0) return(character(0))
  gids <- str_split(enr$geneID[idx[1]], "/")[[1]]
  entrez_map$SYMBOL[match(gids, entrez_map$ENTREZID)]
}

# only store targets here
db_targets <- list(
  KEGG     = kegg_targets,
  Reactome = reactome_targets,
  GO       = go_targets
)

# loop over contrasts
contrasts <- c("Pre_vs_Young","Post_vs_Pre_MA")

for(ctr in contrasts){
  val_col  <- paste0("log2FC_",   ctr)
  pval_col <- paste0("adj_",      ctr)
  
  for(db_name in names(db_targets)){
    # get the right enrichResult object by name
    enrich_res <- get(paste0(tolower(db_name), "_", ctr))
    targets    <- db_targets[[db_name]]
    
    df_all <- bind_rows(lapply(names(targets), function(pathway){
      syms <- get_syms(enrich_res, targets[pathway])
      gene_stats %>% 
        filter(Gene_Symbol %in% syms, .data[[pval_col]] <= 0.10) %>%
        transmute(
          Gene_Symbol,
          value   = .data[[val_col]],
          star    = case_when(
            .data[[pval_col]] < 0.001 ~ "***",
            .data[[pval_col]] < 0.01  ~ "**",
            .data[[pval_col]] < 0.05  ~ "*",
            TRUE                       ~ ""
          ),
          Pathway = pathway
        )
    }))
    if(nrow(df_all)==0) next
    
    df_all <- df_all %>% 
      group_by(Pathway) %>%
      mutate(Gene_Symbol = fct_reorder(Gene_Symbol, value)) %>%
      ungroup()
    
    print(
      ggplot(df_all, aes(x = Gene_Symbol, y = value)) +
        geom_col() +
        geom_text(aes(label = star), hjust = -0.1) +
        coord_flip() +
        facet_wrap(~ Pathway, scales = "free_y", ncol = 2) +
        labs(
          title = paste(db_name, ctr),
          x     = NULL,
          y     = paste0("log₂FC (", ctr, ")")
        ) +
        theme_minimal()
    )
  }
}


# ---- 11. Overlaid Radar Plots for Pathway Clusters ----
library(dplyr)
library(tidyr)
library(fmsb)
library(broom)
library(tibble)
library(ggplot2)
library(stringr)

categories <- list(
  Redox = list(
    Redox_Stress = c("NOX1","NOX2","NOX4","CYBA","CYBB","DUOX1","DUOX2","XDH","AIFM1","AIFM2"),
    Antioxidant_Defenses = c("SOD1","SOD2","SOD3","CAT","GPX1","GPX2","GPX3","GPX4",
                             "PRDX1","PRDX2","PRDX3","PRDX4","PRDX5","PRDX6",
                             "TXN","TXNRD1","TXNRD2","GSR","GCLC","GCLM","HMOX1"),
    Oxidative_Damage_Signalling = c("NFE2L2","KEAP1","NQO1","HMOX1","TXNIP","PARP1","ATM",
                                    "ATR","TP53","NFKB1","RELA")
  ),
  Synthesis = list(
    Protein_Synthesis_Signalling = c("MTOR","RPTOR","RICTOR","AKT1","AKT2","TSC1","TSC2",
                                     "RHEB","PRKAA1","PRKAA2","RPS6KB1","RPS6KB2",
                                     "EIF4EBP1","EIF4E"),
    Ribosome_Biogenesis = c(paste0("RPL", c("P0","P1","P2","3","4","7","11")),
                            paste0("RPS", c("3","6","8","19")),
                            "MRPL12","NOP14","NOP56","UTP15","FBL","PES1","BOP1","WDR12"),
    ER_Stress = c("HSPA5","HSP90B1","CANX","CALR","PDIA1","PDIA3","PDIA4",
                  "PDIA6","DDIT3","GADD34","ATF6","EIF2AK3","ERN1","ERN2",
                  "SEC61A1","DERL1"),
    UPR_and_Chaperones = c("XBP1","XBP1s","ATF4","DDIT3","HSPA5","HSP90AA1","HSP90AB1",
                           "DNAJA1","DNAJB1","DNAJB9","DNAJC3","BAG3","GRP94")
  ),
  Proteolysis = list(
    Calpains = c("CAPN1","CAPN2","CAPN3","CAPN5","CAPN10","CAPNS1","CAST"),
    Caspases = paste0("CASP",  c(1,3,4,5,6,7,8,9,10,12)),
    Autophagy = c("ULK1","BECN1","ATG3","ATG5","ATG7","ATG12","ATG16L1",
                  "MAP1LC3A","MAP1LC3B","SQSTM1"),
    Ubiquitin_Proteasome = c(paste0("PSMA",1:7), paste0("PSMB",1:8),
                             paste0("PSMC",1:6), paste0("PSMD",1:14),
                             "UBB","UBC","UBA1","UBE2D2","UBE2L3","UBE3A")
  ),
  Mito = list(
    Fusion = c("MFN1","MFN2","OPA1","OPA3"),
    Fission = c("DNM1L","FIS1","MFF","MIEF1","MIEF2"),
    Mitophagy = c("PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","BCL2L13"),
    Biogenesis = c("PPARGC1A","PPARGC1B","NRF1","NRF2","TFAM","TFB2M")
  )
)

groups <- c("Young","Pre_MA","Post_MA")
cols   <- c("#1b9e77","#d95f02","#7570b3")
all_abundances <- list()

for (set_name in names(categories)) {
  cats <- categories[[set_name]]
  vars <- names(cats)
  
  # a) Mean abundances
  abund_mat <- sapply(groups, function(g) {
    sub <- long_df %>%
      filter(
        (g=="Young"   & Group=="Young") |
          (g=="Pre_MA"  & Group=="MA"    & Timepoint=="Pre") |
          (g=="Post_MA" & Group=="MA"    & Timepoint=="Post")
      )
    sapply(cats, function(genes) mean(sub$Abundance[sub$Gene_Symbol %in% genes], na.rm=TRUE))
  })
  
  # b) Normalize to Young=1
  abund_norm <- abund_mat / abund_mat[,"Young"]
  max_val    <- ceiling(max(abund_norm, na.rm=TRUE) * 1.1)
  
  # c) Compute module‑level p‑values
  module_p <- tibble(
    Module           = vars,
    p_Pre_vs_Young   = sapply(vars, function(m) {
      dfm <- long_df %>% filter(
        ((Group=="MA" & Timepoint=="Pre") | Group=="Young") &
          Gene_Symbol %in% cats[[m]]
      )
      if (n_distinct(dfm$Group)==2) tidy(t.test(Abundance ~ Group, data=dfm))$p.value else NA
    }),
    p_Post_vs_Young  = sapply(vars, function(m) {
      dfm <- long_df %>% filter(
        ((Group=="MA" & Timepoint=="Post") | Group=="Young") &
          Gene_Symbol %in% cats[[m]]
      )
      if (n_distinct(dfm$Group)==2) tidy(t.test(Abundance ~ Group, data=dfm))$p.value else NA
    }),
    p_Post_vs_Pre_MA = sapply(vars, function(m) {
      dfm <- long_df %>% filter(
        (Group=="MA" & Timepoint %in% c("Pre","Post")) &
          Gene_Symbol %in% cats[[m]]
      )
      if (n_distinct(dfm$Timepoint)==2) tidy(t.test(Abundance ~ Timepoint, data=dfm))$p.value else NA
    })
  )
  
  # d) Print p‑values table
  message("\n=== ", set_name, " module p‑values ===")
  print(module_p)
  
  # e) Store normalized abundances for later
  all_abundances[[set_name]] <-
    as.data.frame(abund_norm, stringsAsFactors = FALSE) %>%
    rownames_to_column("Module") %>%
    pivot_longer(-Module, names_to="Group", values_to="Abundance") %>%
    mutate(Pathway = set_name)
  
  # f) Build radar_df and plot individual radar
  radar_df <- as.data.frame(rbind(
    max = rep(max_val, length(vars)),
    min = rep(0,       length(vars)),
    t(abund_norm)
  ))
  colnames(radar_df) <- vars
  rownames(radar_df) <- c("max","min", groups)
  
  fmsb::radarchart(
    radar_df,
    axistype    = 1,
    seg         = max_val,
    caxislabels = seq(0, max_val, by=1),
    pcol        = cols,
    pfcol       = NA,
    plwd        = 2,
    plty        = 1,
    vlcex       = 0.8,
    title       = paste0(set_name, " (Norm_to_Young)")
  )
  legend(
    "bottomleft",
    legend = groups,
    bty    = "n",
    pch    = 20,
    col    = cols,
    pt.cex = 2,
    cex    = 0.8
  )
}

# ---- Combined Pathway Radar Plot ----
abundances_all_pathways <- bind_rows(all_abundances)
summary_abund <- abundances_all_pathways %>%
  group_by(Pathway, Group) %>%
  summarize(Abundance = mean(Abundance, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=Pathway, values_from=Abundance)

max_comb <- ceiling(max(summary_abund[-1], na.rm=TRUE) * 1.1)
radar_comb <- as.data.frame(rbind(
  max = rep(max_comb, ncol(summary_abund) - 1),
  min = rep(0,       ncol(summary_abund) - 1),
  summary_abund[ , -1]
))
colnames(radar_comb) <- names(summary_abund)[-1]
rownames(radar_comb) <- c("max", "min", summary_abund$Group)

fmsb::radarchart(
  radar_comb,
  axistype    = 1,
  seg         = max_comb,
  caxislabels = seq(0, max_comb, by=1),
  pcol        = cols,
  pfcol       = NA,
  plwd        = 2,
  plty        = 1,
  vlcex       = 0.8,
  title       = "Combined Pathway Abundances (Norm_to_Young)"
)
legend(
  "bottomleft",
  legend = summary_abund$Group,
  bty    = "n",
  pch    = 20,
  col    = cols,
  pt.cex = 2,
  cex    = 0.8
)

# ---- 12. Bar Plots with Correct Group Ordering ----

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Define group order and colors
group_order <- c("Young", "Pre_MA", "Post_MA")
cols <- c("Young" = "#1b9e77", "Pre_MA" = "#d95f02", "Post_MA" = "#7570b3")

# Reconstruct the categories list (must match your earlier definition)
categories <- list(
  Redox = list(
    Redox_Stress = c("NOX1","NOX2","NOX4","CYBA","CYBB","DUOX1","DUOX2","XDH","AIFM1","AIFM2"),
    Antioxidant_Defenses = c("SOD1","SOD2","SOD3","CAT","GPX1","GPX2","GPX3","GPX4",
                             "PRDX1","PRDX2","PRDX3","PRDX4","PRDX5","PRDX6",
                             "TXN","TXNRD1","TXNRD2","GSR","GCLC","GCLM","HMOX1"),
    Oxidative_Damage_Signalling = c("NFE2L2","KEAP1","NQO1","HMOX1","TXNIP","PARP1","ATM",
                                    "ATR","TP53","NFKB1","RELA")
  ),
  Synthesis = list(
    Protein_Synthesis_Signalling = c("MTOR","RPTOR","RICTOR","AKT1","AKT2","TSC1","TSC2",
                                     "RHEB","PRKAA1","PRKAA2","RPS6KB1","RPS6KB2",
                                     "EIF4EBP1","EIF4E"),
    Ribosome_Biogenesis = c(paste0("RPL", c("P0","P1","P2","3","4","7","11")),
                            paste0("RPS", c("3","6","8","19")),
                            "MRPL12","NOP14","NOP56","UTP15","FBL","PES1","BOP1","WDR12"),
    ER_Stress = c("HSPA5","HSP90B1","CANX","CALR","PDIA1","PDIA3","PDIA4",
                  "PDIA6","DDIT3","GADD34","ATF6","EIF2AK3","ERN1","ERN2",
                  "SEC61A1","DERL1"),
    UPR_and_Chaperones = c("XBP1","XBP1s","ATF4","DDIT3","HSPA5","HSP90AA1","HSP90AB1",
                           "DNAJA1","DNAJB1","DNAJB9","DNAJC3","BAG3","GRP94")
  ),
  Proteolysis = list(
    Calpains = c("CAPN1","CAPN2","CAPN3","CAPN5","CAPN10","CAPNS1","CAST"),
    Caspases = paste0("CASP",  c(1,3,4,5,6,7,8,9,10,12)),
    Autophagy = c("ULK1","BECN1","ATG3","ATG5","ATG7","ATG12","ATG16L1",
                  "MAP1LC3A","MAP1LC3B","SQSTM1"),
    Ubiquitin_Proteasome = c(paste0("PSMA",1:7), paste0("PSMB",1:8),
                             paste0("PSMC",1:6), paste0("PSMD",1:14),
                             "UBB","UBC","UBA1","UBE2D2","UBE2L3","UBE3A")
  ),
  Mito = list(
    Fusion = c("MFN1","MFN2","OPA1","OPA3"),
    Fission = c("DNM1L","FIS1","MFF","MIEF1","MIEF2"),
    Mitophagy = c("PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","BCL2L13"),
    Biogenesis = c("PPARGC1A","PPARGC1B","NRF1","NRF2","TFAM","TFB2M")
  )
)

# Prepare module p-values for annotation
all_module_p <- list()
for (set_name in names(categories)) {
  cats <- categories[[set_name]]
  vars <- names(cats)
  module_p <- tibble(
    Module = vars,
    p_Pre_vs_Young = sapply(vars, function(m) {
      dfm <- long_df %>% filter(
        ((Group=="MA" & Timepoint=="Pre") | Group=="Young") &
          Gene_Symbol %in% cats[[m]]
      )
      if (n_distinct(dfm$Group)==2) t.test(Abundance ~ Group, data=dfm)$p.value else NA
    }),
    p_Post_vs_Young = sapply(vars, function(m) {
      dfm <- long_df %>% filter(
        ((Group=="MA" & Timepoint=="Post") | Group=="Young") &
          Gene_Symbol %in% cats[[m]]
      )
      if (n_distinct(dfm$Group)==2) t.test(Abundance ~ Group, data=dfm)$p.value else NA
    })
  )
  all_module_p[[set_name]] <- module_p
}

# Plot each pathway set
for (set_name in names(all_abundances)) {
  plot_df <- all_abundances[[set_name]] %>%
    mutate(Group = factor(Group, levels = group_order))
  
  # Add module p-values for annotation
  plot_df <- plot_df %>%
    left_join(all_module_p[[set_name]], by = "Module")
  
  # Create significance labels
  sig_labels <- plot_df %>%
    mutate(
      label = case_when(
        Group == "Pre_MA" & p_Pre_vs_Young < 0.001 ~ "***",
        Group == "Pre_MA" & p_Pre_vs_Young < 0.01 ~ "**",
        Group == "Pre_MA" & p_Pre_vs_Young < 0.05 ~ "*",
        Group == "Post_MA" & p_Post_vs_Young < 0.001 ~ "***",
        Group == "Post_MA" & p_Post_vs_Young < 0.01 ~ "**",
        Group == "Post_MA" & p_Post_vs_Young < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Bar plot
  p <- ggplot(sig_labels, aes(x = Module, y = Abundance, fill = Group)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = label),
              position = position_dodge(0.8),
              vjust = -0.5,
              size = 4) +
    scale_fill_manual(values = cols) +
    labs(
      title = paste(set_name, "Pathway Abundances"),
      subtitle = "Normalized to Young Group",
      y = "Normalized Abundance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  print(p)
}

# ---- Combined Pathway Bar Plot ----
combined_df <- bind_rows(all_abundances) %>%
  mutate(Group = factor(Group, levels = group_order))

ggplot(combined_df, aes(x = Module, y = Abundance, fill = Group)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Combined Pathway Analysis",
    y = "Normalized Abundance"
  ) +
  facet_wrap(~ Pathway, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
