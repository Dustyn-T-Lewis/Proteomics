# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)

# 1. Load MyoF and SARCO Data from Excel
file_path <- "13_MASTER_SEER proteome analysis (4-4-23).xlsx"
myof <- read_excel(file_path, sheet = "MASTER MyoF data")
sarco <- read_excel(file_path, sheet = "MASTER_SARCO_data", skip = 2)  # skipping rows if needed

# 2. Rename first column for consistency
colnames(myof)[1] <- "Gene_Protein"
colnames(sarco)[1] <- "Gene_Protein"

# 3. Add a column to indicate the fraction for each dataset
myof <- myof %>% mutate(Fraction = "MyoF")
sarco <- sarco %>% mutate(Fraction = "SARCO")

# 4. Select relevant columns (adjust as necessary - here we select columns containing sample data and identifiers)
myof_clean <- myof %>% select(Gene_Protein, Fraction, contains("PRE"), contains("POST"), contains("PPS"))
sarco_clean <- sarco %>% select(Gene_Protein, Fraction, contains("EAA"), contains("PPS"))

# 5. Combine the datasets
combined <- bind_rows(myof_clean, sarco_clean)

# 6. Reshape the combined data into long format:
#    All sample columns become rows, while Gene_Protein and Fraction remain as identifiers.
long_df <- combined %>%
  pivot_longer(cols = -c(Gene_Protein, Fraction),
               names_to = "Sample",
               values_to = "Abundance") %>%
  mutate(
    # Define timepoint based on column names: samples with "PRE" or "POST"
    Timepoint = case_when(
      grepl("PRE", Sample) ~ "Pre",
      grepl("POST", Sample) ~ "Post",
      TRUE ~ NA_character_
    ),
    # Define group based on sample name (example: "PPS" or "EAA" samples as Young, the rest as MA)
    Group = case_when(
      grepl("PPS|EAA", Sample) ~ "Young",
      TRUE ~ "MA"
    )
  )

# 7. Extract the Gene Symbol for PANTHER: take the text before the first underscore in Gene_Protein.
#    This is one of the supported IDs.
long_df <- long_df %>%
  mutate(
    Gene_Symbol = sub("_.*", "", Gene_Protein)
  )

# 8. Display the structure of the working combined dataframe with the new Gene_Symbol column
str(long_df)
head(long_df)

# 9. PANTHER Batch ID Search supported IDs:
#    https://www.pantherdb.org/help/PANTHER_Batch_ID_Search.jsp


      # Now, Compare Young vs Middle-Aged (Pre) using independent t-tests (unpaired).
      #	Compare Middle-Aged Pre vs Post Training using dependent t-tests (paired).
      
      # Significance thresholds:
      # p < 0.01 = significant (especially for protein-level changes)
      # p < 0.05 = reportable (training phenotypes)
      # p < 0.10 = numerical trend due to small sample size
      
      #	Bioinformatics (PANTHER): Proteins were split by direction (up/downregulated) 
      # and analyzed for overrepresentation of GO-Slim biological processes using Fisher’s 
      # test + Bonferroni correction.

# 10. AGE COMPARISON: Independent t-test (Young vs MA) using Pre samples
age_results <- long_df %>%
  filter(Timepoint == "Pre") %>%
  group_by(Gene_Symbol, Fraction) %>%
  filter(!is.na(Abundance)) %>%
  summarize(
    t_test = list(t.test(Abundance ~ Group)),
    .groups = "drop"
  ) %>%
  mutate(
    stats = map(t_test, tidy)
  ) %>%
  unnest(cols = stats) %>%
  select(Gene_Symbol, Fraction, p.value, estimate1, estimate2, statistic, conf.low, conf.high) %>%
  rename(
    Mean_MA = estimate1,
    Mean_Young = estimate2,
    Age_tstat = statistic,
    Age_pval = p.value,
    Age_CI_low = conf.low,
    Age_CI_high = conf.high
  ) %>%
  mutate(
    Age_Significance = case_when(
      Age_pval < 0.01 ~ "Significant",
      Age_pval < 0.05 ~ "Approaching",
      Age_pval < 0.10 ~ "Numerical",
      TRUE ~ "NS"
    )
  )

# 11. TRAINING COMPARISON: Paired t-test (Pre vs Post) in MA group only
training_data <- long_df %>%
  filter(Group == "MA", Timepoint %in% c("Pre", "Post")) %>%
  pivot_wider(names_from = Timepoint, values_from = Abundance)

training_results <- training_data %>%
  group_by(Gene_Symbol, Fraction) %>%
  filter(!is.na(Pre), !is.na(Post)) %>%
  summarize(
    t_test = list(t.test(Pre, Post, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    stats = map(t_test, tidy)
  ) %>%
  unnest(cols = stats) %>%
  select(Gene_Symbol, Fraction, p.value, estimate, statistic, conf.low, conf.high) %>%
  rename(
    Training_Diff = estimate,
    Training_tstat = statistic,
    Training_pval = p.value,
    Training_CI_low = conf.low,
    Training_CI_high = conf.high
  ) %>%
  mutate(
    Training_Significance = case_when(
      Training_pval < 0.01 ~ "Significant",
      Training_pval < 0.05 ~ "Approaching",
      Training_pval < 0.10 ~ "Numerical",
      TRUE ~ "NS"
    )
  )

# 12. Merge both result sets into your working dataframe
final_df <- full_join(age_results, training_results, by = c("Gene_Symbol", "Fraction"))
