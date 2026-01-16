###########################  PREPROCESS THE DATA  ########################### 

# Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(scales)
library(tibble)
library(purrr)

# Read data
df <- read_csv("StringDB_results/ALL_HPO.csv") %>%
  select(where(~ !all(is.na(.)))) 

# Load the updated classification table
categories <- read_csv("HPO_terms_classified.csv") %>%
  select(-term)

# Merge based on 'description'
df <- df %>%  
  left_join(categories, by = "description")

# Exclusions (as provided)
exclude <- c(
  "Abnormal facial shape","Abnormal facial skeleton morphology","Abnormal fetal morphology",
  "Abnormal palate morphology","Abnormal skull morphology","Abnormal soft palate morphology",
  "Abnormality of facial musculature","Abnormality of facial soft tissue",
  "Abnormality of head or neck","Abnormality of prenatal development or birth",
  "Abnormality of skull size","Abnormality of the forehead","Abnormality of the head",
  "Abnormality of the skull base","Aplasia cutis congenita","Aplasia of the fingers",
  "Aplasia/Hypoplasia affecting the eye","Aplasia/Hypoplasia affecting the fundus",
  "Aplasia/Hypoplasia involving bones of the feet",
  "Aplasia/Hypoplasia involving the central nervous system","Aplasia/Hypoplasia of fingers",
  "Aplasia/Hypoplasia of the cerebrum","Aplasia/Hypoplasia of the corpus callosum",
  "Aplasia/Hypoplasia of the mandible","Aplasia/Hypoplasia of the skin",
  "Aplasia/hypoplasia affecting bones of the axial skeleton",
  "Aplasia/hypoplasia involving bones of the extremities",
  "Aplasia/hypoplasia involving bones of the hand",
  "Aplasia/hypoplasia involving bones of the lower limbs",
  "Aplasia/hypoplasia involving bones of the upper limbs","Aplasia/hypoplasia involving the skeleton",
  "Aplasia/hypoplasia of the extremities","Birth weight","Broad forehead","Budd-Chiari syndrome",
  "Cleft palate","Cleft soft palate","Congenital abnormal hair pattern",
  "Congenital malformation of the great arteries","Coronal craniosynostosis","Craniosynostosis",
  "Decreased head circumference","Developmental cataract","Developmental glaucoma",
  "Erythroid hypoplasia","Fetal anomaly","Genetic disorder","Global developmental delay",
  "Head tremor","High palate","\"High, narrow palate\"","Hydrops fetalis",
  "Hypoplasia of teeth","Hypoplasia of the maxilla","Lambdoidal craniosynostosis","Microcephaly",
  "Narrow palate","Neurodevelopmental abnormality","Neurodevelopmental delay",
  "Nonimmune hydrops fetalis","Optic nerve hypoplasia","Oral cleft","Pain in head and neck region",
  "Proportionate short stature","Pure red cell aplasia","Radial artery aplasia","Renal agenesis",
  "Renal hypoplasia/aplasia","Sagittal craniosynostosis","Severe short stature","Short stature",
  "Small forehead","Age","Test result","Sign or symptom","Disease","Sign or symptom",
  "Population measurement","Self reported educational attainment","Parental longevity",
  "Traffic air pollution measurement","Educational attainment","Hypospadias","Narrow mouth",
  "Brachycephaly","Microcephaly","Turricephaly","Abnormality of the face",
  "Abnormal mandible morphology","Retrognathia","Abnormality of the forehead","Pointed chin",
  "Hypertelorism","Triangular face","Abnormality of the maxilla","Broad forehead","Long philtrum",
  "Abnormal location of ears","Abnormal earlobe morphology","Low-set ears","Tooth malposition",
  "Split hand","Triphalangeal thumb","Slender finger","Broad finger","Small for gestational age",
  "Split foot","Metatarsus adductus","Underdeveloped supraorbital ridges","Crumpled ear",
  "Duplication of thumb phalanx","Partial duplication of thumb phalanx",
  "Partial duplication of the phalanx of hand","Duplication of phalanx of hallux",
  "Abnormality of the phalanges of the toes","Midface retrusion","Short digit","Flat face",
  "Spina bifida occulta", "Short palm", "Short 5th finger"
)

# Remove excluded rows
df <- df %>% 
  filter(!description %in% exclude)



###########   GENERATE CHEMICAL CONTRIBUTIONS STACKED BAR CHART   ###########

# Get top 10 categories by absolute count
hpo_top_categories_abs <- df %>%
  count(class, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  pull(class)

# Filter to top 10 categories 
hpo_filtered_abs <- df %>%
  filter(class %in% hpo_top_categories_abs)

# Count (class, ChemicalName)
hpo_stacked_data <- hpo_filtered_abs %>%
  count(class, ChemicalName)

# Reorder categories by total count
hpo_category_order_abs <- hpo_stacked_data %>%
  group_by(class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class)

# Factor with desired order
hpo_stacked_data$class <- factor(
  hpo_stacked_data$class,
  levels = hpo_category_order_abs
)

# Plot stacked bar chart
hpo_chemical_contributions_stacked_bar <- ggplot(hpo_stacked_data, aes(
  x = n,
  y = class,
  fill = ChemicalName
)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Chemical Contributions in Top 10 HPO Categories",
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

print(hpo_chemical_contributions_stacked_bar)

# Save plot
ggsave(
  "HPO_analysis/chemical_contributions_stacked_bar.png",
  plot = hpo_chemical_contributions_stacked_bar,
  width = 10, height = 6, dpi = 300
)



#######   GENERATE RELATIVE CHEMICAL CONTRIBUTIONS STACKED BAR CHART  #######

# Top 10 categories by absolute count (reuse logic for consistency)
hpo_top_categories_rel <- df %>%
  count(class, name = "total") %>%
  slice_max(total, n = 10) %>%
  pull(class)

# Filter to those categories
hpo_filtered_rel <- df %>%
  filter(class %in% hpo_top_categories_rel)

# Count (class, ChemicalName) pairs
hpo_contributions_rel <- hpo_filtered_rel %>%
  count(class, ChemicalName, name = "n")

# Normalize within each category
hpo_normalized_rel <- hpo_contributions_rel %>%
  group_by(class) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

# Order categories by total occurrences (same order as absolute totals)
hpo_category_order_rel <- hpo_contributions_rel %>%
  group_by(class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class)

hpo_normalized_rel$class <- factor(
  hpo_normalized_rel$class,
  levels = hpo_category_order_rel
)

# Plot normalized stacked bar (vertical)
hpo_rel_plot <- ggplot(hpo_normalized_rel, aes(x = class, y = share, fill = ChemicalName)) +
  geom_col(width = 0.7) +
  labs(
    title = "Relative Chemical Contributions in Top 10 HPO Categories",
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

print(hpo_rel_plot)

# Save
ggsave("HPO_analysis/normalized_chemical_contributions.png",
       plot = hpo_rel_plot, width = 12, height = 6, dpi = 300)



############   BLOOD CATEGORY: OCCURRENCE-BASED QUANTIFICATION   ############

# Define HPO groups
blood_chemistry_biomarkers <- c(
  "Protein measurement","Abnormal circulating metabolite concentration","Blood protein measurement",
  "Hematological measurement","Lipid measurement","Lipid or lipoprotein measurement",
  "Triglyceride measurement","Abnormal circulating lipid concentration","Abnormal enzyme/coenzyme activity",
  "Lipoprotein measurement","Abnormal LDL cholesterol concentration",
  "Abnormal circulating carboxylic acid concentration","Abnormal circulating cholesterol concentration",
  "Abnormal circulating nitrogen compound concentration","Apolipoprotein A 1 measurement",
  "Decreased circulating GABA concentration","Elevated sweat chloride","Liver enzyme measurement",
  "Serum alanine aminotransferase measurement","Serum urea measurement",
  "Abnormal blood cation concentration","Abnormal blood ion concentration",
  "Abnormal blood potassium concentration","Abnormal circulating albumin concentration",
  "Abnormal erythrocyte sedimentation rate","Abnormal glucose homeostasis","Cytokine measurement",
  "Glucose measurement","Hypercholesterolemia","Hyperlipoproteinemia",
  "Increased LDL cholesterol concentration","Low density lipoprotein cholesterol measurement",
  "Metabolite measurement","Serum gamma-glutamyl transferase measurement","Transferrin saturation measurement",
  "Abnormal blood monovalent inorganic cation concentration","Abnormal circulating bilirubin concentration",
  "Abnormal circulating creatine kinase concentration","Abnormal circulating dicarboxylic acid concentration",
  "Abnormal circulating ferritin concentration","Abnormal circulating non-proteinogenic amino acid concentration",
  "Amino acid measurement","Antibody measurement","Apolipoprotein B measurement",
  "Aspartate aminotransferase measurement",
  "Aspartate aminotransferase to alanine aminotransferase ratio",
  "Decreased plasma total carnitine","Diabetes mellitus biomarker","Elevated alpha-fetoprotein",
  "Elevated circulating creatine kinase concentration","Endothelial growth factor measurement",
  "Glucose intolerance","High density lipoprotein cholesterol measurement","Hyperammonemia",
  "Hyperlipidemia","Hypoalbuminemia","Iron biomarker measurement",
  "Leucine carboxyl methyltransferase 1 measurement","Liver disease biomarker",
  "Serum albumin measurement","Serum iron measurement","Sphingolipid measurement",
  "Total cholesterol measurement","Urate measurement","Increased circulating ferritin concentration",
  "Elevated erythrocyte sedimentation rate"
)

erythrocytes_rbc <- c(
  "Anemia","Abnormal erythrocyte morphology","Hemoglobin measurement","Anemia of inadequate production",
  "Complete blood cell count","Reticulocyte count","Reticulocyte measurement","Macrocytic anemia",
  "Mean corpuscular volume","Erythrocyte count","Erythrocyte indices","Erythrocyte measurement",
  "Hypochromic anemia","Persistence of hemoglobin F","Red blood cell distribution width",
  "Macrocytic dyserythropoietic anemia","Poikilocytosis","Reticulocytopenia","Methemoglobinemia",
  "Polycythemia","Reticulocytosis","Abnormal number of erythroid precursors",
  "Elevated red cell adenosine deaminase level","Extramedullary hematopoiesis","Heinz bodies",
  "Abnormal hemoglobin","Normochromic anemia","Increased mean corpuscular volume",
  "Mean corpuscular hemoglobin","Mean reticulocyte volume","Abnormal mean corpuscular volume",
  "Red blood cell density measurement"
)

## Percentage of Blood HPO occurrences in all HPO occurrences ##

blood_total_count <- df %>%
  filter(class == "Blood") %>%
  nrow()

blood_in_all_pct <- blood_total_count / nrow(df) * 100

print(sprintf(
  "%.1f%% of all HPO occurrences are classified as Blood (%d out of %d).",
  blood_in_all_pct, blood_total_count, nrow(df)
))

## Percentage of Blood HPO occurrences that are biomarkers ##

biomarkers_count <- df %>%
  filter(class == "Blood") %>%
  filter(grepl(paste(blood_chemistry_biomarkers, collapse = "|"),
               description, ignore.case = TRUE)) %>%
  nrow()

biomarkers_pct <- biomarkers_count / blood_total_count * 100

print(sprintf(
  "%.1f%% of Blood HPO occurrences are biomarkers (%d out of %d).",
  biomarkers_pct, biomarkers_count, blood_total_count
))

## Percentage of Blood HPO occurrences that are erythrocyte/RBC-related ##

rbc_count <- df %>%
  filter(class == "Blood") %>%
  filter(grepl(paste(erythrocytes_rbc, collapse = "|"),
               description, ignore.case = TRUE)) %>%
  nrow()

rbc_pct <- rbc_count / blood_total_count * 100

print(sprintf(
  "%.1f%% of Blood HPO occurrences are erythrocyte/RBC-related (%d out of %d).",
  rbc_pct, rbc_count, blood_total_count
))

## Percentage of Musculoskeletal HPO occurrences in all HPO occurrences ##

musculoskeletal_total_count <- df %>%
  filter(class == "Musculoskeletal") %>%
  nrow()

musculoskeletal_in_all_pct <- musculoskeletal_total_count / nrow(df) * 100

print(sprintf(
  "%.1f%% of all HPO occurrences are classified as Musculoskeletal (%d out of %d).",
  musculoskeletal_in_all_pct, musculoskeletal_total_count, nrow(df)
))

## Percentage of Cardiovascular HPO occurrences in all HPO occurrences ##

cardiovascular_total_count <- df %>%
  filter(class == "Cardiovascular") %>%
  nrow()

cardiovascular_in_all_pct <- cardiovascular_total_count / nrow(df) * 100

print(sprintf(
  "%.1f%% of all HPO occurrences are classified as Cardiovascular (%d out of %d).",
  cardiovascular_in_all_pct, blood_total_count, nrow(df)
))


########################   CNS GROUPING (HPO: BRAIN + MENTAL + SELECTED EYE)   ########################
# NOTE:
# "CNS" stands for Central Nervous System.
# For HPO phenotypes, we group:
#   1) all "Brain" phenotypes
#   2) all "Mental" phenotypes
#   3) ONLY selected "Eye" phenotypes that are plausibly CNS-related (vision pathways + eye movement control),
#      e.g., diplopia/strabismus/abnormal eye movements/abnormal binocular vision/abnormal conjugate movement.
# We do NOT include general/anterior-segment eye phenotypes (e.g., cornea/lens/sclera/conjunctiva) because many
# are structural ocular findings not directly reflecting CNS function.

# Selected Eye phenotypes to include in CNS (from your Eye list)
eye_cns_selected <- c(
  "Diplopia",
  "Abnormality of binocular vision",
  "Abnormality of eye movement",
  "Abnormal conjugate eye movement",
  "Strabismus",
  "Saccadic smooth pursuit"
)

# Create CNS-labeled version of df (DO NOT overwrite df)
df_cns <- df %>%
  mutate(
    class = if_else(
      class %in% c("Brain", "Mental") |
        (class == "Eye" & description %in% eye_cns_selected),
      "CNS",
      class
    )
  )

###########   GENERATE CHEMICAL CONTRIBUTIONS STACKED BAR CHART (CNS VERSION)   ###########

# Get top 10 categories by absolute count
hpo_top_categories_abs_cns <- df_cns %>%
  count(class, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  pull(class)

# Filter to top 10 categories
hpo_filtered_abs_cns <- df_cns %>%
  filter(class %in% hpo_top_categories_abs_cns)

# Count (class, ChemicalName)
hpo_stacked_data_cns <- hpo_filtered_abs_cns %>%
  count(class, ChemicalName)

# Reorder categories by total count
hpo_category_order_abs_cns <- hpo_stacked_data_cns %>%
  group_by(class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class)

# Factor with desired order
hpo_stacked_data_cns$class <- factor(
  hpo_stacked_data_cns$class,
  levels = hpo_category_order_abs_cns
)

# Plot stacked bar chart
hpo_chemical_contributions_stacked_bar_cns <- ggplot(hpo_stacked_data_cns, aes(
  x = n,
  y = class,
  fill = ChemicalName
)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Chemical Contributions in Top 10 HPO Categories (CNS = Brain + Mental + selected Eye)",
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

print(hpo_chemical_contributions_stacked_bar_cns)

# Save plot
ggsave(
  "HPO_analysis/chemical_contributions_stacked_bar_CNS.png",
  plot = hpo_chemical_contributions_stacked_bar_cns,
  width = 10, height = 6, dpi = 300
)



#######   GENERATE RELATIVE CHEMICAL CONTRIBUTIONS STACKED BAR CHART (CNS VERSION)  #######

# Top 10 categories by absolute count
hpo_top_categories_rel_cns <- df_cns %>%
  count(class, name = "total") %>%
  slice_max(total, n = 10) %>%
  pull(class)

# Filter to those categories
hpo_filtered_rel_cns <- df_cns %>%
  filter(class %in% hpo_top_categories_rel_cns)

# Count (class, ChemicalName) pairs
hpo_contributions_rel_cns <- hpo_filtered_rel_cns %>%
  count(class, ChemicalName, name = "n")

# Normalize within each category
hpo_normalized_rel_cns <- hpo_contributions_rel_cns %>%
  group_by(class) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

# Order categories by total occurrences
hpo_category_order_rel_cns <- hpo_contributions_rel_cns %>%
  group_by(class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class)

hpo_normalized_rel_cns$class <- factor(
  hpo_normalized_rel_cns$class,
  levels = hpo_category_order_rel_cns
)

# Plot normalized stacked bar (vertical)
hpo_rel_plot_cns <- ggplot(hpo_normalized_rel_cns, aes(x = class, y = share, fill = ChemicalName)) +
  geom_col(width = 0.7) +
  labs(
    title = "Relative Chemical Contributions in Top 10 HPO Categories (CNS = Brain + Mental + selected Eye)",
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

print(hpo_rel_plot_cns)

# Save
ggsave(
  "HPO_analysis/normalized_chemical_contributions_CNS.png",
  plot = hpo_rel_plot_cns,
  width = 12, height = 6, dpi = 300
)



############   CNS CATEGORY: OCCURRENCE-BASED QUANTIFICATION   ############

# Count all CNS occurrences
cns_total_count <- df_cns %>%
  filter(class == "CNS") %>%
  nrow()

cns_in_all_pct <- cns_total_count / nrow(df_cns) * 100

print(sprintf(
  "%.1f%% of all HPO occurrences are classified as CNS (%d out of %d).",
  cns_in_all_pct, cns_total_count, nrow(df_cns)
))

# Within CNS: breakdown by source class (Brain vs Mental vs Eye-selected)
cns_breakdown <- df_cns %>%
  mutate(
    cns_source = case_when(
      class == "CNS" & (description %in% eye_cns_selected) ~ "Eye (CNS-selected)",
      class == "CNS" & (description %in% unique(df$description[df$class == "Eye"])) ~ "Eye (CNS-selected)",
      TRUE ~ NA_character_
    )
  )

# More robust: compute sources directly from original df
cns_sources <- df %>%
  mutate(
    cns_source = case_when(
      class == "Brain" ~ "Brain",
      class == "Mental" ~ "Mental",
      class == "Eye" & description %in% eye_cns_selected ~ "Eye (CNS-selected)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cns_source))

print(cns_sources %>% count(cns_source, sort = TRUE))

# Percent of CNS occurrences that come from each source
print(
  cns_sources %>%
    count(cns_source, name = "n") %>%
    mutate(pct = n / sum(n) * 100) %>%
    arrange(desc(n))
)

########################   CNS SOURCES PIE CHART (WITH n + pct LABELS)   ########################

cns_pie_data <- cns_sources %>%
  count(cns_source, name = "n") %>%
  mutate(
    pct = n / sum(n) * 100,
    label = paste0(cns_source, "\n", "n = ", n, " | ", percent(pct / 100, accuracy = 0.1))
  ) %>%
  arrange(desc(n))

cns_sources_pie <- ggplot(cns_pie_data, aes(x = "", y = n, fill = cns_source)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4,
            lineheight = 0.95) +
  labs(
    title = "CNS Occurrences by Source (Brain vs Mental vs Eye-selected)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.2),
    legend.title = element_text(size = 12)
  ) +
  guides(fill = guide_legend(title = "CNS Source"))

print(cns_sources_pie)

ggsave(
  "HPO_analysis/cns_sources_pie.png",
  plot = cns_sources_pie,
  width = 8, height = 6, dpi = 300
)

# File 1: Count of HPO descriptions (after filtering)
hpo_description_counts <- df %>%
  count(description, sort = TRUE)

write_csv(
  hpo_description_counts,
  "HPO_analysis/hpo_description_counts.csv"
)

# File 2: Unique inputGenes and preferredNames per HPO description
unique_by_hpo <- df %>%
  select(description, inputGenes, preferredNames) %>%
  mutate(
    inputGenes = if_else(is.na(inputGenes), "", inputGenes),
    preferredNames = if_else(is.na(preferredNames), "", preferredNames)
  ) %>%
  group_by(description) %>%
  summarise(
    inputGenes_unique = {
      vals <- unlist(str_split(inputGenes, ","))
      vals <- str_trim(vals)
      vals <- vals[!is.na(vals) & vals != ""]
      paste(sort(unique(vals)), collapse = ", ")
    },
    preferredNames_unique = {
      vals <- unlist(str_split(preferredNames, ","))
      vals <- str_trim(vals)
      vals <- vals[!is.na(vals) & vals != ""]
      paste(sort(unique(vals)), collapse = ", ")
    },
    .groups = "drop"
  )

write_csv(
  unique_by_hpo,
  "HPO_analysis/hpo_unique_inputGenes_and_preferredNames.csv"
)
