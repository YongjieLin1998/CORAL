# CORAL
CORAL(clone oriented reconstruction of attractors in the landscape) is a project for heritbale cell state identification and gene expression dynamics inference, mainly based on the lineage tracing strategy. 

## CORAL-base

CORAL-base is the original version of CORAL, designed for one timepoint lineage tracing dataset analysis. By integrating lineage information, it provides an alternative way to define the boundary of a really 'heritable' cell state, rather than a 'cluster'. Then the 'CORAL'-like cell state tree enables a systematic analysis of 'memory genes' participating in many biology process. 

Installation
You can install the development version of CORAL from GitHub with:

R

if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("your-github-username/CORAL")
Please replace your-github-username/CORAL with the actual repository URL.

Example Workflow
This guide demonstrates the complete CORAL analysis workflow, from an existing Seurat object to final visualizations.

1. Prerequisites
This tutorial assumes you have already loaded your single-cell data into a Seurat object named seurat_obj. This object must contain:

Normalized expression data (e.g., processed via Seurat::NormalizeData()).

A metadata column (e.g., true_clone_id) that contains the unique lineage barcode or clone identifier for each cell.

First, load the necessary libraries.

R

# Load required libraries
library(CORAL)
library(Seurat)
library(ComplexHeatmap)
library(ggplot2)

# This guide assumes `seurat_obj` is already in your environment.
# Let's inspect its metadata to confirm the required columns are present.
# The output should show columns like 'true_clone_id' and 'cell_type'.
head(seurat_obj@meta.data)
2. Run the Core CORAL Analysis
The run_coral_ground_truth_analysis() function is the core of the workflow. It calculates clone distances, defines CORAL states, and identifies heritable genes. All results are conveniently stored within the Seurat object.

R

# Run the main CORAL analysis pipeline with 4 states
seurat_obj <- run_coral_ground_truth_analysis(
  seurat_obj = seurat_obj,
  true_barcode_col = "true_clone_id", # Specify the clone ID column
  num_states = 4,                      # Group clones into 4 states
  permutation_repeats = 100,           # Recommended repeats for robust results
  n_cores = 4                          # Number of cores for parallel processing
)

# The results are stored in the 'misc' slot of the Seurat object
results_4_states <- seurat_obj@misc$CORAL_ground_truth_analysis
3. Iterating and Modifying the Number of States
A key part of the analysis is exploring the ideal number of states to describe your clonal hierarchy. You can easily test a different number of states by re-running the main analysis function with a new value for num_states. This will overwrite the previous results.

R

# Let's say we want to explore what 6 states look like instead of 4
message("Re-running analysis with 6 states...")
seurat_obj <- run_coral_ground_truth_analysis(
  seurat_obj = seurat_obj,
  true_barcode_col = "true_clone_id",
  num_states = 6, # <-- Changed value
  permutation_repeats = 100,
  n_cores = 4
)

# Now, all downstream visualizations will use the new 6-state grouping.
# For example, let's regenerate the MDS plot with updated states:
p_mds_6_states <- visualize_clone_mds(seurat_obj, color_by = "coral_state")
print(p_mds_6_states)
4. Run Advanced Gene Fluctuation Analysis
The analyze_gene_fluctuation() function provides deeper insight into the dynamics of heritable genes. It should be run after you have settled on a final number of states.

R

# Run the gene fluctuation mode analysis on the 6-state result
seurat_obj <- analyze_gene_fluctuation(seurat_obj, n_cores = 4)

# The results are updated in the heritable_genes_df with a new 'omega_area' column
print("Heritable genes dataframe with omega_area:")
head(seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df)
5. Visualize the Final Results
CORAL provides a suite of plotting functions to explore the final analysis output.

Clone Energy Distance Heatmap
This heatmap displays the transcriptional similarity between all clones. The clustering reveals clone groupings, which are now partitioned into 6 CORAL states.

R

# Visualize the clone-clone energy distance matrix with 6 states
ht <- visualize_clone_distance_heatmap(seurat_obj)
ComplexHeatmap::draw(ht)
Heritable Gene Distribution
This plot is a key diagnostic tool that compares the observed distribution of Omega-squared values (heritability effect size) against the null distribution from permutations.

R

# Compare the observed vs. null distribution of Omega-squared
p_dist <- plot_heritable_gene_distribution(seurat_obj)
print(p_dist)
Confusion Matrix: CORAL States vs. Cell Types
This heatmap shows the correspondence between our 6 lineage-defined CORAL states and any other pre-existing cell annotation, such as cell type.

R

# Create a confusion matrix heatmap
# Assumes a 'cell_type' column exists in your metadata
p_confusion <- plot_state_celltype_confusion(
  seurat_obj,
  celltype_col = "cell_type"
)
ComplexHeatmap::draw(p_confusion)
Gene Fluctuation Mode Plot
This plot visualizes the relationship between a gene's heritability effect size (Omega-squared) and its fluctuation pattern (omega_area), revealing genes with distinct modes of inheritance.

R

# Visualize the gene fluctuation mode
p_fluctuation <- plot_gene_fluctuation_mode(
  seurat_obj,
  genes_to_highlight = c("GATA1", "PAX5", "CD34") # Highlight genes of interest
)
print(p_fluctuation)
(Note: The figures/ directory is a suggestion. Please generate these plots and place them in your repository, updating the paths accordingly.)

Citation
If you use CORAL in your research, please cite: Lin, Y., Chen, X., Wu, L., Zhou, Y., & Lin, Y. (2025). Widespread transcriptional memory shapes heritable states and functional heterogeneity in cancer and stem cells. bioRxiv, 2025-08.
The original codes and example datasets of the [initial manuscript](https://www.biorxiv.org/content/10.1101/2025.08.21.671653v1.full) are available in https://drive.google.com/drive/folders/1-cNiSKZFyVSs9Mndq87AcRXfaGweLesj?usp=sharing. 


