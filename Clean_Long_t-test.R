#### Proteomic Data Processing & T-Test Pipeline ####

# ---- Load Required Libraries ----
library(readxl)     # Read Excel files
library(dplyr)      # Data wrangling
library(stringr)    # String operations
library(tidyr)      # Reshaping data
library(httr)       # HTTP requests
library(jsonlite)   # JSON parsing
library(purrr)      # Functional programming

# ---- 1. Load & Clean Data ----
sarco <- read_excel("13_MASTER_SEER proteome analysis (4-4-23).xlsx", sheet = "MASTER_SARCO_data") %>%
  select(where(~ any(!is.na(.)))) %>%
  filter(if_any(everything(), ~ !is.na(.))) %>%
  mutate(
    Gene_Symbol = str_extract(Gene_protein, "^[^_;]+") %>% str_remove(";?GN=+$"),
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
    j <- content(r, as = "parsed", type = "application/json")
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
  select(-Gene_Symbol_API)

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
  select(-matches("^T-test_|^DELTA_|^CORREL_|^Young_|^detection_")) %>%   # Safely exclude unwanted statistical columns
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
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

#### 6. Statistical Comparison ####

# ---- 6a. Unpaired T-Tests: Young vs MA Pre ----
results1 <- intersect(
  unique(long_df$Gene_Symbol[long_df$Group == "Young"]),
  unique(long_df$Gene_Symbol[long_df$Timepoint == "Pre"])
) %>%
  map_dfr(function(gene) {
    y <- long_df %>% filter(Group == "Young", Gene_Symbol == gene) %>% pull(Abundance)
    pre <- long_df %>% filter(Group == "MA", Timepoint == "Pre", Gene_Symbol == gene) %>% pull(Abundance)
    if (length(y) > 1 && length(pre) > 1) {
      t <- t.test(y, pre)
      tibble(
        Gene_Symbol = gene,
        Young_Mean = mean(y), Pre_MA_Mean = mean(pre),
        Change = mean(pre) - mean(y),
        p_value = t$p.value
      )
    }
  })

# ---- 6b. T-Tests: MA Pre vs MA Post ----
results2 <- intersect(
  unique(long_df$Gene_Symbol[long_df$Timepoint == "Pre"]),
  unique(long_df$Gene_Symbol[long_df$Timepoint == "Post"])
) %>%
  map_dfr(function(gene) {
    pre <- long_df %>% filter(Group == "MA", Timepoint == "Pre", Gene_Symbol == gene) %>% pull(Abundance)
    post <- long_df %>% filter(Group == "MA", Timepoint == "Post", Gene_Symbol == gene) %>% pull(Abundance)
    if (length(pre) > 1 && length(post) > 1) {
      t <- t.test(pre, post)
      tibble(
        Gene_Symbol = gene,
        Pre_MA_Mean = mean(pre), Post_MA_Mean = mean(post),
        Change = mean(post) - mean(pre),
        p_value = t$p.value
      )
    }
  })

# ---- Output Results ----
print(results1)
print(results2)

#### 7. Bioinformatics via Panther ####

# FILTER by p < 0.01 for highly significant protein results
sig1 <- results1 %>% filter(p_value < 0.01)  # Young vs MA Pre
sig2 <- results2 %>% filter(p_value < 0.01)  # MA Pre vs Post

# FILTER by 0.01 ≤ p < 0.05 for "approaching" significance
sig3 <- results1 %>% filter(p_value >= 0.01 & p_value < 0.05)
sig4 <- results2 %>% filter(p_value >= 0.01 & p_value < 0.05)


# SPLIT SIGNIFICANT into UP and DOWN regulated genes
results1_up_significant   <- sig1 %>% filter(Change > 0) %>% pull(Gene_Symbol) %>% unique()
results1_down_significant <- sig1 %>% filter(Change < 0) %>% pull(Gene_Symbol) %>% unique()

results2_up_significant   <- sig2 %>% filter(Change > 0) %>% pull(Gene_Symbol) %>% unique()
results2_down_significant <- sig2 %>% filter(Change < 0) %>% pull(Gene_Symbol) %>% unique()

# SPLIT APPROACHING SIGNIFICANCE into UP and DOWN
results1_up_approaching   <- sig3 %>% filter(Change > 0) %>% pull(Gene_Symbol) %>% unique()
results1_down_approaching <- sig3 %>% filter(Change < 0) %>% pull(Gene_Symbol) %>% unique()

results2_up_approaching   <- sig4 %>% filter(Change > 0) %>% pull(Gene_Symbol) %>% unique()
results2_down_approaching <- sig4 %>% filter(Change < 0) %>% pull(Gene_Symbol) %>% unique()

# WRITE to text files for PANTHER (SIGNIFICANT RESULTS)
write.table(results1_up_significant,   "Y_vs_Pre_UP_significant.txt",   row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(results1_down_significant, "Y_vs_Pre_DOWN_significant.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(results2_up_significant,   "Pre_vs_Post_UP_significant.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(results2_down_significant, "Pre_vs_Post_DOWN_significant.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)

# WRITE to text files for PANTHER (APPROACHING RESULTS)
write.table(results1_up_approaching,   "Y_vs_Pre_UP_approaching.txt",   row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(results1_down_approaching, "Y_vs_Pre_DOWN_approaching.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(results2_up_approaching,   "Pre_vs_Post_UP_approaching.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(results2_down_approaching, "Pre_vs_Post_DOWN_approaching.txt",row.names = FALSE, col.names = FALSE, quote = FALSE)


# Ready for subsequent PANTHER analysis

#### Gene Enrichment of gene sets ####
young_df <- results1 %>%
  select(Gene_Symbol, Abundance = Young_Mean)

pre_df <- results1 %>%
  select(Gene_Symbol, Abundance = Pre_MA_Mean)

post_df <- results2 %>%
  select(Gene_Symbol, Abundance = Post_MA_Mean)

# Write to CSV files
write_tsv(young_df, "young_mean_abundance.txt")
write_tsv(pre_df, "pre_ma_mean_abundance.txt")
write_tsv(post_df, "post_ma_mean_abundance.txt")
