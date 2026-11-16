# Comparative 16S rRNA Microbiome Analysis

## Project Overview
This project presents a comparative analysis of 16S rRNA microbiome communities from three different environments:

- Gut
- Marine sediment
- Soil

The aim is to compare microbial diversity, community composition, and differences in microbial abundance across these environments.

## Dataset
The dataset contains 18 samples:
- 6 gut samples
- 6 marine sediment samples
- 6 soil samples

## Analysis Methods
The analysis was performed using R and includes:

1. Data cleaning and preparation
2. Alpha diversity analysis
3. Beta diversity analysis
4. Principal Coordinates Analysis (PCoA)
5. PERMANOVA statistical testing
6. Differential abundance analysis using DESeq2
7. Taxonomic composition analysis

## Repository Contents

| File | Description |
|---|---|
| `Complete_Microbiome_Analysis.R` | Complete R analysis script |
| `Final_Microbiome_Analysis.Rmd` | R Markdown source for the report |
| `Final_Microbiome_Analysis_Academic_Review.pdf` | Final academic report |
| `.png` files | Diversity, taxonomic composition, PCoA, and biomarker figures |
| `.csv` files | Statistical results, summaries, and analysis tables |

## Software and Packages
The analysis was conducted using R. Packages include:

- phyloseq
- vegan
- DESeq2
- ggplot2
- dplyr
- tidyr
- readxl
- rmarkdown

## How to Use
1. Download or clone this repository.
2. Open the project in RStudio.
3. Place the required input dataset in the appropriate `Data` folder.
4. Check that the file paths in the R script match your project folders.
5. Install the required R packages.
6. Run `Complete_Microbiome_Analysis.R` in RStudio.

**Note:** The original input dataset is not included in this repository. It must be supplied separately to reproduce the analysis.

## Report
The final academic report is available as:

`Final_Microbiome_Analysis_Academic_Review.pdf`

---

This repository was prepared for academic review and demonstrates a comparative 16S rRNA microbiome analysis of gut, marine sediment, and soil communities.
