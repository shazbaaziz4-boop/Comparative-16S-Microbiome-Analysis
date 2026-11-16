
# ============================================================
# COMPLETE MICROBIOME ANALYSIS
# ============================================================
# Project: Comparative analysis of Gut, Marine and Soil samples
# Total samples: 18 (6 samples per group)
# Main methods:
# - Data preparation and filtering
# - Alpha diversity
# - Beta diversity
# - Taxonomic composition
# - DESeq2 differential abundance analysis
#
# ============================================================
# SECTION 1: PROJECT SETUP
# ============================================================

# Set random seed for reproducibility
set.seed(123)

# Check project working directory
getwd()

# Required R packages
required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "phyloseq",
  "vegan",
  "DESeq2",
  "tibble"
)

# Check which packages are available
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace,
          logical(1), quietly = TRUE)
]

# Display missing packages
print(missing_packages)

# Create output folders if they do not exist
dir.create("Results", showWarnings = FALSE)
dir.create("Report", showWarnings = FALSE)
dir.create("Report/Figures", recursive = TRUE,
           showWarnings = FALSE)
dir.create("Report/Tables", recursive = TRUE,
           showWarnings = FALSE)

# ============================================================
# END OF PROJECT SETUP
# ============================================================


# ============================================================
# SECTION 2: LOAD EXISTING PHYLOSEQ OBJECTS
# ============================================================

# Load the cleaned phyloseq object
ps_clean <- readRDS("Results/cleaned_phyloseq.rds")

# Load the rarefied phyloseq object
ps_rarefied <- readRDS("Results/rarefied_phyloseq.rds")

# Check the cleaned dataset
print(ps_clean)

# Check the rarefied dataset
print(ps_rarefied)

# Check the number of samples in each group
print(table(sample_data(ps_clean)$Group))

# Check sequencing depth after rarefaction
print(range(sample_sums(ps_rarefied)))

# ============================================================
# END OF SECTION 2
# ============================================================



# ============================================================
# SECTION 3: ALPHA DIVERSITY ANALYSIS
# ============================================================

# Load packages required by this script
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

# Calculate alpha diversity using the rarefied dataset
alpha_results <- estimate_richness(
  ps_rarefied,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)

# Add sample IDs
alpha_results$SampleID <- rownames(alpha_results)

# Prepare sample metadata
alpha_metadata <- data.frame(
  SampleID = sample_names(ps_rarefied),
  Group = as.character(sample_data(ps_rarefied)$Group)
)

# Combine diversity results with sample groups
alpha_results <- left_join(
  alpha_results,
  alpha_metadata,
  by = "SampleID"
)

# Inspect the results
print(head(alpha_results))
print(table(alpha_results$Group))

# Save alpha diversity results
write.csv(
  alpha_results,
  "Results/Alpha_Diversity_Results.csv",
  row.names = FALSE
)

# Calculate group summaries
alpha_group_summary <- alpha_results %>%
  group_by(Group) %>%
  summarise(
    Observed_Median = median(Observed),
    Observed_IQR = IQR(Observed),
    Chao1_Median = median(Chao1),
    Chao1_IQR = IQR(Chao1),
    Shannon_Median = median(Shannon),
    Shannon_IQR = IQR(Shannon),
    Simpson_Median = median(Simpson),
    Simpson_IQR = IQR(Simpson),
    .groups = "drop"
  )

print(alpha_group_summary)

write.csv(
  alpha_group_summary,
  "Report/Tables/Alpha_Group_Summary.csv",
  row.names = FALSE
)

# Perform Kruskal-Wallis tests
alpha_measures <- c(
  "Observed", "Chao1", "Shannon", "Simpson"
)

alpha_kw_results <- lapply(alpha_measures, function(measure) {
  test <- kruskal.test(
    alpha_results[[measure]] ~ alpha_results$Group
  )

  data.frame(
    Measure = measure,
    Statistic = unname(test$statistic),
    P_value = test$p.value
  )
})

alpha_kw_results <- bind_rows(alpha_kw_results)

print(alpha_kw_results)

write.csv(
  alpha_kw_results,
  "Report/Tables/Alpha_Kruskal_Wallis.csv",
  row.names = FALSE
)


# Perform pairwise Wilcoxon tests with BH correction
alpha_pairwise_results <- lapply(alpha_measures, function(measure) {
  
  test <- pairwise.wilcox.test(
    alpha_results[[measure]],
    alpha_results$Group,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  p_matrix <- test$p.value
  
  valid_pairs <- which(!is.na(p_matrix), arr.ind = TRUE)
  
  data.frame(
    Measure = measure,
    Comparison = paste(
      rownames(p_matrix)[valid_pairs[, 1]],
      "vs",
      colnames(p_matrix)[valid_pairs[, 2]]
    ),
    Adjusted_P_value = p_matrix[valid_pairs]
  )
})

alpha_pairwise_results <- bind_rows(alpha_pairwise_results)

print(alpha_pairwise_results)

write.csv(
  alpha_pairwise_results,
  "Report/Tables/Alpha_Pairwise_Summary.csv",
  row.names = FALSE
)


# Generate and save four alpha diversity boxplots
for (measure in alpha_measures) {

  p <- ggplot(
    alpha_results,
    aes(x = Group, y = .data[[measure]], fill = Group)
  ) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.12, size = 2) +
    theme_bw() +
    labs(
      title = paste(measure, "Diversity"),
      x = "Sample Group",
      y = measure
    ) +
    theme(legend.position = "none")

  ggsave(
    filename = paste0("Report/Figures/", measure, "_Diversity.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# END OF SECTION 3
# ============================================================



# ============================================================
# SECTION 4: BETA DIVERSITY ANALYSIS
# ============================================================

library(phyloseq)
library(vegan)
library(ggplot2)

# Convert cleaned counts to relative abundance
ps_tss <- transform_sample_counts(
  ps_clean,
  function(x) x / sum(x)
)

# Extract OTU table
otu_tss <- as(otu_table(ps_tss), "matrix")

# Ensure samples are rows
if (taxa_are_rows(ps_tss)) {
  otu_tss <- t(otu_tss)
}

# Prepare metadata
metadata_beta <- data.frame(
  Group = as.factor(sample_data(ps_tss)$Group)
)

rownames(metadata_beta) <- sample_names(ps_tss)

# ------------------------------------------------------------
# 4.1 Bray-Curtis dissimilarity and PCoA
# ------------------------------------------------------------

bray_dist <- vegdist(
  otu_tss,
  method = "bray"
)

pcoa_result <- cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE
)

bray_variance <- round(
  pcoa_result$eig[1:2] /
    sum(pcoa_result$eig[pcoa_result$eig > 0]) * 100,
  2
)

bray_plot_data <- data.frame(
  SampleID = rownames(pcoa_result$points),
  PCoA1 = pcoa_result$points[, 1],
  PCoA2 = pcoa_result$points[, 2]
)

bray_plot_data$Group <- metadata_beta[
  bray_plot_data$SampleID, "Group"
]

pcoa_plot <- ggplot(
  bray_plot_data,
  aes(PCoA1, PCoA2, color = Group)
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCoA - Bray-Curtis Dissimilarity",
    x = paste0("PCoA1 (", bray_variance[1], "%)"),
    y = paste0("PCoA2 (", bray_variance[2], "%)")
  )

ggsave(
  "Report/Figures/PCoA_BrayCurtis.png",
  pcoa_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4.2 Jaccard dissimilarity and PCoA
# ------------------------------------------------------------

jaccard_dist <- vegdist(
  otu_tss,
  method = "jaccard",
  binary = TRUE
)

jaccard_pcoa <- cmdscale(
  jaccard_dist,
  k = 2,
  eig = TRUE
)

jaccard_variance <- round(
  jaccard_pcoa$eig[1:2] /
    sum(jaccard_pcoa$eig[jaccard_pcoa$eig > 0]) * 100,
  2
)

jaccard_plot_data <- data.frame(
  SampleID = rownames(jaccard_pcoa$points),
  PCoA1 = jaccard_pcoa$points[, 1],
  PCoA2 = jaccard_pcoa$points[, 2]
)

jaccard_plot_data$Group <- metadata_beta[
  jaccard_plot_data$SampleID, "Group"
]

jaccard_plot <- ggplot(
  jaccard_plot_data,
  aes(PCoA1, PCoA2, color = Group)
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCoA - Jaccard Dissimilarity",
    x = paste0("PCoA1 (", jaccard_variance[1], "%)"),
    y = paste0("PCoA2 (", jaccard_variance[2], "%)")
  )

ggsave(
  "Report/Figures/PCoA_Jaccard.png",
  jaccard_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4.3 PERMANOVA tests
# ------------------------------------------------------------

permanova_bray <- adonis2(
  bray_dist ~ Group,
  data = metadata_beta,
  permutations = 999
)

permanova_jaccard <- adonis2(
  jaccard_dist ~ Group,
  data = metadata_beta,
  permutations = 999
)

print(permanova_bray)
print(permanova_jaccard)

write.csv(
  as.data.frame(permanova_bray),
  "Report/Tables/PERMANOVA_BrayCurtis.csv"
)

write.csv(
  as.data.frame(permanova_jaccard),
  "Report/Tables/PERMANOVA_Jaccard.csv"
)

# ------------------------------------------------------------
# 4.4 Test homogeneity of group dispersion
# ------------------------------------------------------------

bray_dispersion <- betadisper(
  bray_dist,
  metadata_beta$Group
)

jaccard_dispersion <- betadisper(
  jaccard_dist,
  metadata_beta$Group
)

print(permutest(bray_dispersion, permutations = 999))
print(permutest(jaccard_dispersion, permutations = 999))

# ============================================================
# END OF SECTION 4
# ============================================================



# ============================================================
# SECTION 4: BETA DIVERSITY ANALYSIS
# ============================================================

library(phyloseq)
library(vegan)
library(ggplot2)

# Convert cleaned counts to relative abundance
ps_tss <- transform_sample_counts(
  ps_clean,
  function(x) x / sum(x)
)

# Extract OTU table
otu_tss <- as(otu_table(ps_tss), "matrix")

# Ensure samples are rows
if (taxa_are_rows(ps_tss)) {
  otu_tss <- t(otu_tss)
}

# Prepare metadata
metadata_beta <- data.frame(
  Group = as.factor(sample_data(ps_tss)$Group)
)

rownames(metadata_beta) <- sample_names(ps_tss)

# ------------------------------------------------------------
# 4.1 Bray-Curtis dissimilarity and PCoA
# ------------------------------------------------------------

bray_dist <- vegdist(
  otu_tss,
  method = "bray"
)

pcoa_result <- cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE
)

bray_variance <- round(
  pcoa_result$eig[1:2] /
    sum(pcoa_result$eig[pcoa_result$eig > 0]) * 100,
  2
)

bray_plot_data <- data.frame(
  SampleID = rownames(pcoa_result$points),
  PCoA1 = pcoa_result$points[, 1],
  PCoA2 = pcoa_result$points[, 2]
)

bray_plot_data$Group <- metadata_beta[
  bray_plot_data$SampleID, "Group"
]

pcoa_plot <- ggplot(
  bray_plot_data,
  aes(PCoA1, PCoA2, color = Group)
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCoA - Bray-Curtis Dissimilarity",
    x = paste0("PCoA1 (", bray_variance[1], "%)"),
    y = paste0("PCoA2 (", bray_variance[2], "%)")
  )

ggsave(
  "Report/Figures/PCoA_BrayCurtis.png",
  pcoa_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4.2 Jaccard dissimilarity and PCoA
# ------------------------------------------------------------

jaccard_dist <- vegdist(
  otu_tss,
  method = "jaccard",
  binary = TRUE
)

jaccard_pcoa <- cmdscale(
  jaccard_dist,
  k = 2,
  eig = TRUE
)

jaccard_variance <- round(
  jaccard_pcoa$eig[1:2] /
    sum(jaccard_pcoa$eig[jaccard_pcoa$eig > 0]) * 100,
  2
)

jaccard_plot_data <- data.frame(
  SampleID = rownames(jaccard_pcoa$points),
  PCoA1 = jaccard_pcoa$points[, 1],
  PCoA2 = jaccard_pcoa$points[, 2]
)

jaccard_plot_data$Group <- metadata_beta[
  jaccard_plot_data$SampleID, "Group"
]

jaccard_plot <- ggplot(
  jaccard_plot_data,
  aes(PCoA1, PCoA2, color = Group)
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCoA - Jaccard Dissimilarity",
    x = paste0("PCoA1 (", jaccard_variance[1], "%)"),
    y = paste0("PCoA2 (", jaccard_variance[2], "%)")
  )

ggsave(
  "Report/Figures/PCoA_Jaccard.png",
  jaccard_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 4.3 PERMANOVA tests
# ------------------------------------------------------------

permanova_bray <- adonis2(
  bray_dist ~ Group,
  data = metadata_beta,
  permutations = 999
)

permanova_jaccard <- adonis2(
  jaccard_dist ~ Group,
  data = metadata_beta,
  permutations = 999
)

print(permanova_bray)
print(permanova_jaccard)

write.csv(
  as.data.frame(permanova_bray),
  "Report/Tables/PERMANOVA_BrayCurtis.csv"
)

write.csv(
  as.data.frame(permanova_jaccard),
  "Report/Tables/PERMANOVA_Jaccard.csv"
)

# ------------------------------------------------------------
# 4.4 Test homogeneity of group dispersion
# ------------------------------------------------------------

bray_dispersion <- betadisper(
  bray_dist,
  metadata_beta$Group
)

jaccard_dispersion <- betadisper(
  jaccard_dist,
  metadata_beta$Group
)

print(permutest(bray_dispersion, permutations = 999))
print(permutest(jaccard_dispersion, permutations = 999))

# ============================================================
# END OF SECTION 4
# ============================================================


# ============================================================
# SECTION 5: TAXONOMIC COMPOSITION
# ============================================================

# The taxonomic composition figures were generated
# previously using the original relative-abundance tables.

# Define the five taxonomic levels
taxonomic_levels <- c(
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus"
)

# Check that all five figures exist
composition_files <- file.path(
  "Report/Figures",
  paste0(taxonomic_levels, "_Composition.png")
)

composition_check <- data.frame(
  Taxonomic_Level = taxonomic_levels,
  Figure_File = composition_files,
  File_Exists = file.exists(composition_files)
)

print(composition_check)

# Stop if any expected figure is missing
if (!all(composition_check$File_Exists)) {
  stop("One or more taxonomic composition figures are missing.")
}

# Save a record of the figure inventory
write.csv(
  composition_check,
  "Report/Tables/Taxonomic_Figure_Inventory.csv",
  row.names = FALSE
)

# ============================================================
# END OF SECTION 5
# ============================================================



# ============================================================
# SECTION 6: DESeq2 DIFFERENTIAL ABUNDANCE ANALYSIS
# ============================================================

library(phyloseq)
library(DESeq2)
library(dplyr)

# Use the cleaned, unrarefied count data
# Gut is the reference group
sample_data(ps_clean)$Group <- factor(
  sample_data(ps_clean)$Group,
  levels = c("Gut", "Marine", "Soil")
)

# Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(
  ps_clean,
  ~ Group
)

# Remove very low-abundance OTUs
dds <- dds[
  rowSums(counts(dds)) >= 10,
]

# Run DESeq2
dds <- DESeq(dds)

# Inspect model coefficients
print(resultsNames(dds))

# ------------------------------------------------------------
# 6.1 Extract results for all three comparisons
# ------------------------------------------------------------

results_marine_gut <- results(
  dds,
  contrast = c("Group", "Marine", "Gut")
)

results_soil_gut <- results(
  dds,
  contrast = c("Group", "Soil", "Gut")
)

results_marine_soil <- results(
  dds,
  contrast = c("Group", "Marine", "Soil")
)

# ------------------------------------------------------------
# 6.2 Prepare taxonomy information
# ------------------------------------------------------------

taxonomy_df <- as.data.frame(
  tax_table(ps_clean)
)

taxonomy_df$OTU_ID <- rownames(taxonomy_df)

# Function to extract significant OTUs and annotate taxonomy
prepare_deseq_results <- function(
    result_object,
    comparison,
    first_group,
    second_group) {
  
  result_df <- as.data.frame(result_object)
  
  result_df$OTU_ID <- rownames(result_df)
  
  result_df <- result_df %>%
    filter(
      !is.na(padj),
      padj < 0.05
    )
  
  result_df$Direction <- ifelse(
    result_df$log2FoldChange > 0,
    paste("Higher in", first_group),
    paste("Higher in", second_group)
  )
  
  result_df$Comparison <- comparison
  
  result_df <- left_join(
    result_df,
    taxonomy_df,
    by = "OTU_ID"
  )
  
  result_df
}

# Prepare annotated results
significant_marine_gut_annotated <-
  prepare_deseq_results(
    results_marine_gut,
    "Marine vs Gut",
    "Marine",
    "Gut"
  )

significant_soil_gut_annotated <-
  prepare_deseq_results(
    results_soil_gut,
    "Soil vs Gut",
    "Soil",
    "Gut"
  )

significant_marine_soil_annotated <-
  prepare_deseq_results(
    results_marine_soil,
    "Marine vs Soil",
    "Marine",
    "Soil"
  )

# ------------------------------------------------------------
# 6.3 Save significant OTU results
# ------------------------------------------------------------

write.csv(
  significant_marine_gut_annotated,
  "Results/Significant_Marine_vs_Gut_Annotated.csv",
  row.names = FALSE
)

write.csv(
  significant_soil_gut_annotated,
  "Results/Significant_Soil_vs_Gut_Annotated.csv",
  row.names = FALSE
)

write.csv(
  significant_marine_soil_annotated,
  "Results/Significant_Marine_vs_Soil_Annotated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6.4 Summarise significant OTUs
# ------------------------------------------------------------

deseq_summary <- data.frame(
  Comparison = c(
    "Marine vs Gut",
    "Soil vs Gut",
    "Marine vs Soil"
  ),
  
  Significant_OTUs = c(
    nrow(significant_marine_gut_annotated),
    nrow(significant_soil_gut_annotated),
    nrow(significant_marine_soil_annotated)
  ),
  
  Higher_in_First_Group = c(
    sum(significant_marine_gut_annotated$Direction == "Higher in Marine"),
    sum(significant_soil_gut_annotated$Direction == "Higher in Soil"),
    sum(significant_marine_soil_annotated$Direction == "Higher in Marine")
  ),
  
  Higher_in_Second_Group = c(
    sum(significant_marine_gut_annotated$Direction == "Higher in Gut"),
    sum(significant_soil_gut_annotated$Direction == "Higher in Gut"),
    sum(significant_marine_soil_annotated$Direction == "Higher in Soil")
  )
)

print(deseq_summary)

write.csv(
  deseq_summary,
  "Report/Tables/DESeq2_Biomarker_Summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6.5 Summarise significant OTUs by phylum
# ------------------------------------------------------------

all_significant_deseq <- bind_rows(
  significant_marine_gut_annotated,
  significant_soil_gut_annotated,
  significant_marine_soil_annotated
)

phylum_biomarker_summary <- all_significant_deseq %>%
  mutate(
    Phylum = ifelse(
      is.na(Phylum) | Phylum == "",
      "Unclassified",
      as.character(Phylum)
    )
  ) %>%
  group_by(Comparison, Phylum, Direction) %>%
  summarise(
    Number_of_OTUs = n(),
    .groups = "drop"
  )

print(head(phylum_biomarker_summary))

write.csv(
  phylum_biomarker_summary,
  "Report/Tables/DESeq2_Phylum_Summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6.6 Copy annotated results to report tables
# ------------------------------------------------------------

file.copy(
  c(
    "Results/Significant_Marine_vs_Gut_Annotated.csv",
    "Results/Significant_Soil_vs_Gut_Annotated.csv",
    "Results/Significant_Marine_vs_Soil_Annotated.csv"
  ),
  "Report/Tables",
  overwrite = TRUE
)

# ============================================================
# END OF SECTION 6
# ============================================================
