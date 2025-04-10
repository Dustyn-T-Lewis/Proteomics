# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(httr)
library(jsonlite)
library(broom)
library(purrr)

#### Formatting Data ####
# 1. Load MyoF and SARCO Data from Excel
file_path <- "13_MASTER_SEER proteome analysis (4-4-23).xlsx"
sarco <- read_excel(file_path, sheet = "MASTER_SARCO_data")  


# 2. Remove empty rows and columns
sarco <- sarco %>% select(where(~ any(!is.na(.))))
sarco <- sarco %>% filter(if_any(everything(), ~ !is.na(.)))

# 3. Extract Gene Symbol from first part of string before underscore or semicolon
sarco <- sarco %>%
  mutate(
    Gene_Symbol = str_extract(Gene_protein, "^[^_;]+"), # looks at each Gene_protein value and grabs 
                                                        # everything before ";" or"_" 
    Gene_Symbol = str_remove(Gene_Symbol, ";?GN=+$") # clean stray suffix (";GN=") 
  )

#### 4. Obtaining gene symbols for values in 'Gene_protein' that lack gene symbols, but have accession numbers "buried"

# 4.1 Extract UniProt accession numbers (e.g., Q9UHX1, P31994)
sarco <- sarco %>%
  mutate(Accession = str_extract(Gene_protein, "\\b[A-Z][0-9][A-Z0-9]{3}[0-9]\\b"))

# 4.2 Prepare unique accessions for query (only where gene symbol is missing)
sample_cols <- grep("^(EAA|PRE|POST|PPS)", names(sarco), value = TRUE)

accessions_to_query <- sarco %>% 
  mutate(detection_score = rowSums(across(all_of(sample_cols)), na.rm = TRUE)) %>% # Compute detection_score
  group_by(Accession, Gene_Symbol) %>%                                              # Group by Accession & Gene_Symbol
  slice_max(order_by = detection_score, n = 1, with_ties = FALSE) %>%                # Keep highest detection per group
  ungroup() %>%
  filter(is.na(Gene_Symbol) | Gene_Symbol == "") %>%                                  # Filter rows missing Gene_Symbol
  pull(Accession) %>%                                                                # Extract Accession values
  unique() %>%                                                                     # Remove duplicate accessions
  na.omit()

print(accessions_to_query)

# 4.3 Define function to query UniProt API for gene symbols
query_uniprot <- function(accessions) {
  url <- "https://rest.uniprot.org/uniprotkb/search"
  query <- paste(accessions, collapse = " OR ")
  
  response <- GET(url, query = list(
    query = query,
    fields = "accession,gene_primary",
    format = "json",
    size = 500
  ))
  
  if (status_code(response) == 200) {
    results <- content(response, as = "parsed", type = "application/json")
    if (!is.null(results$results)) {
      tibble(
        Accession = sapply(results$results, function(x) x$primaryAccession),
        Gene_Symbol_API = sapply(results$results, function(x) {
          if (!is.null(x$genes) && length(x$genes) > 0 &&
              !is.null(x$genes[[1]]$geneName$value)) {
            x$genes[[1]]$geneName$value
          } else {
            NA_character_
          }
        })
      )
    } else {
      tibble(Accession = character(), Gene_Symbol_API = character())
    }
  } else {
    stop("UniProt query failed with status: ", status_code(response))
  }
}

# 4.4 Batch accessions and query
uniprot_results <- accessions_to_query %>%
  split(ceiling(seq_along(.) / 100)) %>%  # Split into batches of 100
  lapply(query_uniprot) %>%               # Query each batch
  bind_rows() %>%                         # Combine results
  distinct(Accession, .keep_all = TRUE)   # Ensure one entry per Accession

# 4.5 Merge results and intelligently update gene symbols
existing_symbols <- unique(na.omit(sarco$Gene_Symbol))  # Already present primary gene symbols

sarco <- sarco %>%
  left_join(uniprot_results, by = "Accession") %>%   # Merge on Accession
  mutate(
    # If Gene_Symbol is missing, update based on whether Gene_protein indicates an isoform.
    # If Gene_protein contains the word "isoform" (case-insensitive), prefix the UniProt gene symbol with "isoform_".
    Gene_Symbol = case_when(
      (is.na(Gene_Symbol) | Gene_Symbol == "") & str_detect(Gene_protein, regex("isoform", ignore_case = TRUE)) ~ paste0("isoform_", Gene_Symbol_API),
      (is.na(Gene_Symbol) | Gene_Symbol == "") & !str_detect(Gene_protein, regex("isoform", ignore_case = TRUE)) ~ Gene_Symbol_API,
      TRUE ~ Gene_Symbol
    )
  ) %>%
  select(-Gene_Symbol_API)

# 4.6 Remove duplicate rows by retaining the best-detected entry and remove isoform rows if a primary exists

sarco_clean <- sarco %>%
  mutate(
    detection_score = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
    # Create a base symbol by removing any "isoform_" prefix
    base_symbol = str_remove(Gene_Symbol, "^isoform_")
  ) %>%
  # First, for each (Accession, base_symbol), keep only the row with the highest detection score
  group_by(Accession, base_symbol) %>%
  slice_max(order_by = detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Then, within each base gene group, if at least one primary (non-isoform) exists, remove rows labeled as isoform
  group_by(base_symbol) %>%
  mutate(has_primary = any(!str_detect(Gene_Symbol, "^isoform_"))) %>%
  ungroup() %>%
  filter(!(has_primary & str_detect(Gene_Symbol, "^isoform_"))) %>%
  select(-detection_score, -base_symbol, -has_primary)

# sarco_clean now contains unique, high-detection entries.
# Primary entries (without "isoform_") are retained if they exist;
# isoform-labeled rows are removed if a primary exists.


# 5. Reshape the combined data into long format:
#    All sample columns become rows, while Gene_Symbol and Fraction remain as identifiers.
long_df <- sarco_clean %>%
  pivot_longer(
    cols = -c(Gene_protein, Gene_Symbol, Accession), 
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  mutate(
    # Define timepoint if desired
    Timepoint = case_when(
      str_detect(Sample, "PRE") ~ "Pre",
      str_detect(Sample, "POST") ~ "Post",
      TRUE ~ NA_character_
    ),
    # Define group if desired
    Group = case_when(
      str_detect(Sample, "EAA|PPS") ~ "Young",
      TRUE ~ "MA"
    )
  )

#### Panther Analysis ####
# Now, Compare Young vs Middle-Aged (Pre) using independent t-tests (unpaired).
#	Compare Middle-Aged Pre vs Post Training using dependent t-tests (paired).

# Significance thresholds:
  # p < 0.01 = significant (especially for protein-level changes)
  # p < 0.05 = reportable (training phenotypes)
  # p < 0.10 = numerical trend due to small sample size

#	Bioinformatics (PANTHER): Proteins were split by direction (up/downregulated) 
# and analyzed for overrepresentation of GO-Slim biological processes using Fisher’s 
# test + Bonferroni correction.

# PANTHER Batch ID Search supported IDs: https://www.pantherdb.org/help/PANTHER_Batch_ID_Search.jsp
  # Long_df is my processed data with the following columns:
    #   - Gene_protein, Gene_Symbol, Accession (identifiers)
    #   - Sample (e.g., "EAA30", "PRE_S14", "POST_S21", etc.)
    #   - Abundance (numeric measurements)
    #   - Timepoint ("Pre" or "Post")
    #   - Group (e.g., "Young" for EAA/PPS samples, otherwise "MA")
            
  # For training analyses, we assume that for MA-group samples, the "Sample" names are
  # of the form "PRE_<Subject>" and "POST_<Subject>" so that we can extract a subject ID.
  # Here, we extract the subject ID as the part after the underscore in the sample name.
# AGE ANALYSIS: Independent T-Tests (Young vs. MA at Pre) ####
