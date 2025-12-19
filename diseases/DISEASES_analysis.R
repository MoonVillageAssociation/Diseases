
###########################  PREPROCESS THE DATA  ########################### 

# Load the libraries

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(forcats)
library(scales)

# Read data
df <- read_csv("StringDB_results/ALL_DISEASES.csv")
class_map <- read_csv("DISEASES_terms_classified.csv")

df <- df %>%
  select(-any_of("term"))

class_map <- class_map %>%
  mutate(
    class = str_replace_all(class, "\\s*;\\s*", ";")
  )

# Join and coalesce categories and remove genetic disorders
df <- df %>%
  left_join(class_map, by = "description", suffix = c(".db", ".map")) %>%
  mutate(
    class = str_replace_all(class, "\\s*;\\s*", ";")
  ) %>%
  select(-any_of(c("class.map", "class.db"))) %>%
  filter(`class` != "Genetic Disorder") 

# Expand multi-class rows 
df_expanded <- df %>%
  separate_rows(class, sep = ";") %>%
  filter(!is.na(class), class != "") 




###########   GENERATE CHEMICAL CONTRIBUTIONS STACKED BAR CHART   ###########


# Get top 10 categories 
top_categories <- df_expanded %>%
  count(class, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  pull(class)

# Filter to top 10 categories 
filtered_df <- df_expanded %>%
  filter(class %in% top_categories)

# Count (class, ChemicalName)
stacked_data <- filtered_df %>%
  count(class, ChemicalName)

# Reorder categories by total count
class_order <- stacked_data %>%
  group_by(class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class)

# Convert class column to factor with correct order
stacked_data$class <- factor(stacked_data$class, levels = class_order)

# Plot stacked bar chart
chemical_contributions_stacked_bar <- ggplot(stacked_data, aes(
  x = n,
  y = class,
  fill = ChemicalName
)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Chemical Contributions in Top 10 Disease Categories",
    x = "Number of Occurrences",
    y = "class"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  guides(fill = guide_legend(title = "Chemical Name"))

print(chemical_contributions_stacked_bar)

# Save plot
ggsave(
  "DISEASES_analysis/chemical_contributions_stacked_bar.png",
  plot = chemical_contributions_stacked_bar,
  width = 10, height = 6, dpi = 300
)




#######   GENERATE RELATIVE CHEMICAL CONTRIBUTIONS STACKED BAR CHART  #######


# Top 10 categories by absolute count
top_categories_rel <- df_expanded %>%
  count(class, name = "total") %>%
  slice_max(total, n = 10) %>%
  pull(class)

# Filter to those categories
filtered_rel <- df_expanded %>%
  filter(class %in% top_categories_rel)

# Count (class, ChemicalName) pairs
contributions_rel <- filtered_rel %>%
  count(class, ChemicalName, name = "n")

# Normalize within each class
normalized_rel <- contributions_rel %>%
  group_by(class) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

# Order categories by total occurrences (same order as absolute totals)
class_order_rel <- contributions_rel %>%
  group_by(class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class)

normalized_rel$class <- factor(normalized_rel$class, levels = class_order_rel)

# Plot normalized stacked bar (vertical)
p_rel <- ggplot(normalized_rel, aes(x = class, y = share, fill = ChemicalName)) +
  geom_col(width = 0.7) +
  labs(
    title = "Relative Chemical Contributions in Top 10 Disease Categories",
    x = "class",
    y = "Relative Share"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12)
  ) +
  guides(fill = guide_legend(title = "Chemical Name"))

print(p_rel)

# Save
ggsave("DISEASES_analysis/normalized_chemical_contributions.png", plot = p_rel, width = 12, height = 6, dpi = 300)




#######   NEUROLOGICAL DISORDERS IN-DEPTH ANALYSIS  #######


#######   NEUROLOGICAL DISORDERS IN-DEPTH ANALYSIS  #######
# NOTE: df_expanded = occurrence-level (each disease–class assignment is one row)

## Calculate percentage of neurological disorders in all disorders ##

# Count total neurological disorders
neurological_total_count <- df_expanded %>%
  filter(class == "Neurological Disease") %>%
  nrow()

neuro_in_all_pct <- neurological_total_count / nrow(df) * 100

print(sprintf(
  "%.1f%% of all diseases are classified as neurological disorders (%d out of %d).",
  neuro_in_all_pct, neurological_total_count, nrow(df)
))


# Define the list of disease descriptions linked to protein misfolding disorders (PMDs)
protein_formation_descriptions <- c(
  "Prion disease", "Primary cutaneous amyloidosis", 
  "APP-related cerebral amyloid angiopathy", 
  "Alzheimers disease", "Creutzfeldt–Jakob disease", 
  "Lewy body dementia", "Neurodegenerative disease", 
  "ITM2B-related cerebral amyloid angiopathy 2", 
  "CST3-related cerebral amyloid angiopathy", 
  "ITM2B-related cerebral amyloid angiopathy 1", 
  "Parkinsons disease", "Familial visceral amyloidosis", 
  "Amyloidosis", "Dementia", "Vascular dementia",
  "Frontotemporal dementia"
)


## Calculate percentage of neurological disorders related to PMDs ##

# Count how many neurological disorders are related to PMDs
protein_formation_in_neuro_count <- df_expanded %>%
  filter(class == "Neurological Disease") %>%
  filter(grepl(
    paste(protein_formation_descriptions, collapse = "|"),
    description,
    ignore.case = TRUE
  )) %>%
  nrow()

# Calculate the proportion
pmd_in_neuro_pct <- protein_formation_in_neuro_count / neurological_total_count * 100

# Print result
print(sprintf(
  "%.1f%% of neurological disorders occurences are related to protein misfolding disorders (PMDs) (%d out of %d occurences).",
  pmd_in_neuro_pct, protein_formation_in_neuro_count, neurological_total_count
))


## Calculate percentage of all diseases related to PMDs ##

# Count all diseases related to PMDs
protein_formation_total_count <- df_expanded %>%
  filter(grepl(
    paste(protein_formation_descriptions, collapse = "|"),
    description,
    ignore.case = TRUE
  )) %>%
  nrow()

# Calculate the proportion
pmd_in_all_pct <- protein_formation_total_count / nrow(df) * 100

# Print result
print(sprintf(
  "%.1f%% of all disease occurences are related to protein misfolding disorders (PMDs) (%d out of %d occurences).",
  pmd_in_all_pct, protein_formation_total_count, nrow(df)
))




#######   HEMATOLOGIC DISORDERS IN-DEPTH ANALYSIS  #######


# Count all hematologic diseases 
hematologic_total_count <- df_expanded %>%
  filter(class == "Hematologic Disease") %>%
  nrow()

hema_in_all_pct <- hematologic_total_count / nrow(df) * 100

print(sprintf(
  "%.1f%% of all diseases are classified as hematologic disorders (%d out of %d).",
  hema_in_all_pct, hematologic_total_count, nrow(df)
))

# Percentage of hematologic diseases that are anemia 
anemia_count <- df_expanded %>%
  filter(class == "Hematologic Disease") %>%
  filter(grepl("anemia", description, ignore.case = TRUE)) %>%
  nrow()

anemia_pct <- anemia_count / hematologic_total_count * 100

print(sprintf(
  "%.1f%% of hematologic diseases are a type of anemia (%d out of %d).",
  anemia_pct, anemia_count, hematologic_total_count
))


# Percentage of hematologic diseases that are cancers 
hematologic_cancer_terms <- c(
  "Hematologic cancer", "Leukemia", "Lymphatic system cancer",
  "Lymphoma", "non-Hodgkin lymphoma",
  "Mature T-cell and NK-cell lymphoma"
)

hematologic_cancer_count <- df_expanded %>%
  filter(class == "Hematologic Disease") %>%
  filter(grepl(
    paste(hematologic_cancer_terms, collapse = "|"),
    description,
    ignore.case = TRUE
  )) %>%
  nrow()

hematologic_cancer_pct <- hematologic_cancer_count / hematologic_total_count * 100

print(sprintf(
  "%.1f%% of hematologic diseases are cancers (%d out of %d).",
  hematologic_cancer_pct, hematologic_cancer_count, hematologic_total_count
))

# Print total hematologic cancer occurrences (explicitly)
print(sprintf(
  "Total hematologic cancer occurrences (subset of hematologic): %d.",
  hematologic_cancer_count
))


#######   CANCERS IN-DEPTH ANALYSIS  #######


## Total cancer occurrences
cancer_total_count <- df_expanded %>%
  filter(class == "Cancer") %>%
  nrow()

cancer_in_all_pct <- cancer_total_count / nrow(df) * 100

print(sprintf(
  "%.1f%% of all diseases are classified as cancers (%d out of %d).",
  cancer_in_all_pct, cancer_total_count, nrow(df)
))

# Percentage of hematologic cancers among all cancers 
hematologic_cancer_terms <- c(
  "Hematologic cancer",
  "Leukemia",
  "Lymphatic system cancer",
  "Lymphoma",
  "Mature T-cell and NK-cell lymphoma",
  "non-Hodgkin lymphoma"
)

## Hematological cancer occurrences among all cancers
hematologic_cancers_count <- df_expanded %>%
  filter(class == "Cancer") %>%
  filter(description %in% hematologic_cancer_terms) %>%
  nrow()

hematologic_among_all_cancers_pct <- hematologic_cancers_count / cancer_total_count * 100

print(sprintf(
  "%.1f%% of all cancer occurrences are hematological cancers (%d out of %d).",
  hematologic_among_all_cancers_pct, hematologic_cancers_count, cancer_total_count
))

