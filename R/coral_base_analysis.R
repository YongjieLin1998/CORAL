################################################################################
#
#       CORAL Package: Ground-Truth Analysis Module (coral-base)
#
# This file contains a suite of functions for analyzing single-cell datasets
# that have known, ground-truth lineage tracing information.
#
#
################################################################################

#' @import Seurat
#' @import ggplot2
#' @import Matrix
#' @importFrom Matrix rowSums colSums
#' @importFrom stats cmdscale quantile aggregate IQR median na.omit
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
#' @importFrom ggrepel geom_label_repel
#' @importFrom grDevices pdf dev.off
NULL

# A hack to prevent NOTE on check for global variables
utils::globalVariables(c("."))


#' @keywords internal
bigSparseDist <- function(x) {
  x <- as.matrix(x)
  rs2 <- Matrix::rowSums(x^2)
  rs2_matrix <- outer(rs2, rs2, "+")
  x_prod <- tcrossprod(x, 2 * x)
  sqrt(abs(rs2_matrix - x_prod))
}

#' @keywords internal
bigSparseDist_pairwise <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  rs2_x <- Matrix::rowSums(x^2)
  rs2_y <- Matrix::rowSums(y^2)
  rs2_matrix <- outer(rs2_x, rs2_y, "+")
  x_prod <- tcrossprod(x, 2 * y)
  sqrt(abs(rs2_matrix - x_prod))
}


#' @title Run the Complete (Optimized) CORAL-base Ground Truth Analysis
#' @description This is the main wrapper function for the `coral.base` analysis module.
#'   This optimized version avoids creating the N x N cell-cell distance matrix in memory,
#'   making it scalable to large datasets. It analyzes the structure of known lineage
#'   tracing clones, defines CORAL states, and identifies "heritable genes". All
#'   results are stored in the Seurat object's misc slot under `CORAL_ground_truth_analysis`.
#'
#' @param seurat_obj A Seurat object containing true clone IDs in its metadata.
#' @param true_barcode_col A string specifying the column name in the metadata that contains the true lineage barcode IDs.
#' @param clone_size_cutoff An integer. Clones smaller than this size will be excluded from the analysis. Defaults to `2`.
#' @param num_states An integer. The number of CORAL states (i.e., clone clusters) to create. Defaults to `6`.
#' @param hclust_method A string specifying the hierarchical clustering method to use. Defaults to `"ward.D"`. See `?hclust` for more options.
#' @param min_gene_mean A numeric value used to filter out lowly expressed genes in the heritable gene analysis. Defaults to `0.05`.
#' @param permutation_repeats An integer specifying the number of repeats for the permutation test to determine the Omega-squared significance threshold. Defaults to `5`.
#' @param n_cores An integer specifying the number of cores to use for parallel processing. Defaults to `parallel::detectCores() - 1`.
#' @param down_sample A numeric value between 0 and 1 for downsampling clones, primarily for quick testing. Defaults to `1` (no downsampling).
#'
#' @return An updated Seurat object with the analysis results stored in `@misc$CORAL_ground_truth_analysis`.
#'
#' @seealso \code{\link{analyze_gene_fluctuation}}, \code{\link{visualize_clone_distance_heatmap}}, \code{\link{visualize_clone_mds}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. Create a mock Seurat object for demonstration
#' counts <- matrix(rpois(1000 * 200, 0.5), ncol = 200, nrow = 1000)
#' rownames(counts) <- paste0("Gene", 1:1000)
#' colnames(counts) <- paste0("Cell", 1:200)
#' metadata <- data.frame(
#'   row.names = colnames(counts),
#'   # Create 20 clones, each with 10 cells
#'   true_clone_id = rep(paste0("Clone", 1:20), each = 10),
#'   cell_type = sample(c("A", "B"), 200, replace = TRUE)
#' )
#' seurat_obj <- CreateSeuratObject(counts, meta.data = metadata)
#' seurat_obj <- NormalizeData(seurat_obj)
#'
#' # 2. Run the main CORAL ground truth analysis
#' seurat_obj <- run_coral_ground_truth_analysis(
#'   seurat_obj = seurat_obj,
#'   true_barcode_col = "true_clone_id",
#'   num_states = 4,
#'   permutation_repeats = 3, # Lower for quick example
#'   n_cores = 2 # Use 2 cores for the example
#' )
#'
#' # 3. Run the advanced gene fluctuation analysis
#' seurat_obj <- analyze_gene_fluctuation(seurat_obj, n_cores = 2)
#'
#' # 4. Visualize the results
#'
#' # Clone distance heatmap
#' ht <- visualize_clone_distance_heatmap(seurat_obj)
#' ComplexHeatmap::draw(ht)
#'
#' # MDS plot of clones, colored by CORAL state
#' mds_plot <- visualize_clone_mds(seurat_obj, color_by = "coral_state")
#' print(mds_plot)
#'
#' # Distribution of heritable genes
#' dist_plot <- plot_heritable_gene_distribution(seurat_obj)
#' print(dist_plot)
#'
#' # CORAL states on UMAP (requires running UMAP first)
#' # seurat_obj <- FindVariableFeatures(seurat_obj)
#' # seurat_obj <- ScaleData(seurat_obj)
#' # seurat_obj <- RunPCA(seurat_obj)
#' # seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
#' # umap_plot <- visualize_coral_states_umap(seurat_obj)
#' # print(umap_plot)
#'
#' # Confusion matrix of states vs. cell types
#' confusion_heatmap <- plot_state_celltype_confusion(
#'   seurat_obj,
#'   celltype_col = "cell_type"
#' )
#' ComplexHeatmap::draw(confusion_heatmap)
#'
#' # Gene fluctuation mode plot
#' fluctuation_plot <- plot_gene_fluctuation_mode(
#'   seurat_obj,
#'   genes_to_highlight = c("Gene1", "Gene5", "Gene10")
#' )
#' print(fluctuation_plot)
#' }
#' @title Run the Complete (Optimized) CORAL-base Ground Truth Analysis (v4 - Robust)
#' @description This is the complete, robustly fixed version of the main CORAL analysis function.
#' It includes a sanity check for metadata/matrix cell consistency and safeguards against
#' data dimension reduction errors for single-cell clones.
#'
#' @param seurat_obj A Seurat object containing true clone IDs in its metadata.
#' @param true_barcode_col A string specifying the column name in the metadata that contains the true lineage barcode IDs.
#' @param clone_size_cutoff An integer. Clones smaller than this size will be excluded from the analysis. Defaults to `2`.
#' @param num_states An integer. The number of CORAL states (i.e., clone clusters) to create. Defaults to `6`.
#' @param hclust_method A string specifying the hierarchical clustering method to use. Defaults to `"ward.D"`. See `?hclust` for more options.
#' @param min_gene_mean A numeric value used to filter out lowly expressed genes in the heritable gene analysis. Defaults to `0.05`.
#' @param permutation_repeats An integer specifying the number of repeats for the permutation test. Defaults to `5`.
#' @param n_cores An integer specifying the number of cores to use for parallel processing. Defaults to `parallel::detectCores() - 1`.
#' @param down_sample A numeric value between 0 and 1 for downsampling clones. Defaults to `1` (no downsampling).
#'
#' @return An updated Seurat object with the analysis results stored in `@misc$CORAL_ground_truth_analysis`.
#' @export
#' @title Run the Complete (Optimized) CORAL-base Ground Truth Analysis (v4 - Robust)
#' @description This is the complete, robustly fixed version of the main CORAL analysis function.
#' It includes a sanity check for metadata/matrix cell consistency and safeguards against
#' data dimension reduction errors for single-cell clones.
#'
#' @param seurat_obj A Seurat object containing true clone IDs in its metadata.
#' @param true_barcode_col A string specifying the column name in the metadata that contains the true lineage barcode IDs.
#' @param clone_size_cutoff An integer. Clones smaller than this size will be excluded from the analysis. Defaults to `2`.
#' @param num_states An integer. The number of CORAL states (i.e., clone clusters) to create. Defaults to `6`.
#' @param hclust_method A string specifying the hierarchical clustering method to use. Defaults to `"ward.D"`. See `?hclust` for more options.
#' @param min_gene_mean A numeric value used to filter out lowly expressed genes in the heritable gene analysis. Defaults to `0.05`.
#' @param permutation_repeats An integer specifying the number of repeats for the permutation test. Defaults to `5`.
#' @param n_cores An integer specifying the number of cores to use for parallel processing. Defaults to `parallel::detectCores() - 1`.
#' @param down_sample A numeric value between 0 and 1 for downsampling clones. Defaults to `1` (no downsampling).
#'
#' @return An updated Seurat object with the analysis results stored in `@misc$CORAL_ground_truth_analysis`.
#' @export
run_coral_ground_truth_analysis <- function(
    seurat_obj,
    true_barcode_col,
    clone_size_cutoff = 2,
    num_states = 6,
    hclust_method = "ward.D",
    min_gene_mean = 0.05,
    permutation_repeats = 5,
    n_cores = NULL,
    down_sample = 1
) {
  # --- 1. Prepare Barcode Data ---
  message("Step 1/6: Preparing barcode data...")
  if (!true_barcode_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Barcode column '", true_barcode_col, "' not found in Seurat object metadata."))
  }
  seurat_obj@meta.data[[true_barcode_col]] <- as.character(seurat_obj@meta.data[[true_barcode_col]])
  barcode_before_filter <- seurat_obj@meta.data[[true_barcode_col]]
  barcode_before_filter[is.na(barcode_before_filter)] <- "0"
  unique_barcodes <- unique(barcode_before_filter)
  barcode_map <- setNames(seq_along(unique_barcodes), unique_barcodes)
  barcode <- unname(barcode_map[barcode_before_filter])
  barcode <- as.factor(barcode)
  seurat_obj <- AddMetaData(seurat_obj, barcode, col.name = "coral_barcode_numeric")
  barcode_sta_df <- as.data.frame(sort(table(barcode), decreasing = TRUE))
  colnames(barcode_sta_df) <- c("BarcodeID", "Freq")
  unlabeled_numeric_id <- barcode_map["0"]
  if (!is.na(unlabeled_numeric_id)) {
    barcode_sta_df <- barcode_sta_df[barcode_sta_df$BarcodeID != as.character(unlabeled_numeric_id), ]
  }
  if (down_sample < 1 && down_sample > 0) {
    message(paste("Downsampling clones to", down_sample * 100, "%..."))
    barcode_sta_df <- barcode_sta_df[sample(nrow(barcode_sta_df), ceiling(down_sample * nrow(barcode_sta_df))), ]
    barcode_sta_df <- barcode_sta_df[order(barcode_sta_df$Freq, decreasing = TRUE), ]
  }

  # --- 2. OPTIMIZED Clone-Clone Energy Distance Calculation ---
  message("Step 2/6: Calculating clone-clone energy distance (Memory-Optimized)...")
  
  # Initial filtering based on metadata counts
  efficient_clones_initial <- barcode_sta_df$BarcodeID[barcode_sta_df$Freq >= clone_size_cutoff]
  if(length(efficient_clones_initial) < 1) stop("No clones found with size >= clone_size_cutoff based on metadata.")
  
  efficient_clone_index <- seurat_obj$coral_barcode_numeric %in% efficient_clones_initial
  
  if(sum(efficient_clone_index) == 0) {
      stop("No cells remain after filtering for efficient clones. Check for inconsistencies between metadata and assay data.")
  }

  sc_expression <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")[, efficient_clone_index, drop = FALSE]
  barcode_used <- seurat_obj$coral_barcode_numeric[efficient_clone_index]

  # --- [ROBUSTNESS FIX]: Sanity check and re-filter clones ---
  message("  - Verifying clone sizes against the expression matrix...")
  true_freq_in_matrix <- as.data.frame(table(barcode_used))
  colnames(true_freq_in_matrix) <- c("BarcodeID", "TrueFreq")
  
  barcode_sta_df_merged <- merge(barcode_sta_df, true_freq_in_matrix, by = "BarcodeID", all.x = TRUE)
  barcode_sta_df_merged$TrueFreq[is.na(barcode_sta_df_merged$TrueFreq)] <- 0
  
  final_efficient_clones <- barcode_sta_df_merged$BarcodeID[barcode_sta_df_merged$TrueFreq >= clone_size_cutoff]
  
  final_clone_number <- length(final_efficient_clones)
  if(final_clone_number < 2) stop(paste("Fewer than 2 clones remaining after verifying sizes against the expression matrix. Initial clones:", length(efficient_clones_initial), "Final clones:", final_clone_number))
  
  message(paste("  - Initial clones passing cutoff:", length(efficient_clones_initial), "| Final valid clones after verification:", final_clone_number))
  
  final_index_in_barcode_used <- barcode_used %in% final_efficient_clones
  sc_expression <- sc_expression[, final_index_in_barcode_used, drop = FALSE]
  barcode_used <- droplevels(barcode_used[final_index_in_barcode_used]) # droplevels is good practice
  barcode_sta_df_filtered <- subset(barcode_sta_df_merged, BarcodeID %in% final_efficient_clones)
  # --- [END OF FIX] ---

  index_map <- split(seq_along(barcode_used), barcode_used)
  clone_indices <- lapply(as.character(barcode_sta_df_filtered$BarcodeID), function(id) index_map[[id]])
  denominators <- barcode_sta_df_filtered$TrueFreq * (barcode_sta_df_filtered$TrueFreq - 1)
  denominators[denominators <= 0] <- 1

  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  message("  - Calculating distances on-the-fly in parallel...")
  clone_dist_rows <- foreach::foreach(
      i = 1:final_clone_number,
      .packages = "Matrix",
      .export = c("bigSparseDist", "bigSparseDist_pairwise")
  ) %dopar% {
    row <- numeric(final_clone_number)
    sub_matrix_i_T <- t(sc_expression[, clone_indices[[i]], drop = FALSE])

    for (j in i:final_clone_number) {
      if (i == j) {
        if (nrow(sub_matrix_i_T) > 1) {
          row[j] <- sum(bigSparseDist(sub_matrix_i_T)) / denominators[i]
        } else {
          row[j] <- 0
        }
      } else {
        sub_matrix_j_T <- t(sc_expression[, clone_indices[[j]], drop = FALSE])
        row[j] <- mean(bigSparseDist_pairwise(sub_matrix_i_T, sub_matrix_j_T))
      }
    }
    row
  }
  parallel::stopCluster(cl)
  
  clone_dist <- matrix(0, nrow = final_clone_number, ncol = final_clone_number)
    # The result from foreach is a list of rows, so we combine them into a matrix
    if (is.list(clone_dist_rows)) {
        clone_dist <- do.call(rbind, clone_dist_rows)
    } else {
        # If it's already a matrix (e.g., from older foreach versions)
        clone_dist <- clone_dist_rows
    }

  clone_dist[lower.tri(clone_dist)] <- t(clone_dist)[lower.tri(clone_dist)]

  diag_clone <- diag(clone_dist)
  clone_dist_E <- 2 * clone_dist - outer(diag_clone, diag_clone, "+")
  rownames(clone_dist_E) <- as.character(barcode_sta_df_filtered$BarcodeID)
  colnames(clone_dist_E) <- as.character(barcode_sta_df_filtered$BarcodeID)

  # --- 3. Define CORAL States ---
  message("Step 3/6: Defining CORAL states...")
  ht <- ComplexHeatmap::Heatmap(clone_dist_E, row_dend_reorder = TRUE, column_dend_reorder = TRUE,
                                clustering_method_rows = hclust_method, clustering_method_columns = hclust_method,
                                row_split = num_states, column_split = num_states)
  ht_drawn <- ComplexHeatmap::draw(ht)
  lineage_order <- ComplexHeatmap::row_order(ht_drawn)

  barcode_cluster <- rep(0, length(seurat_obj$coral_barcode_numeric))
  names(barcode_cluster) <- colnames(seurat_obj)
  for (i in seq_along(lineage_order)) {
    target_barcodes <- as.character(barcode_sta_df_filtered$BarcodeID[lineage_order[[i]]])
    cells_in_cluster_indices <- which(seurat_obj$coral_barcode_numeric %in% target_barcodes)
    barcode_cluster[cells_in_cluster_indices] <- i
  }
  barcode_cluster <- as.factor(barcode_cluster)
  seurat_obj <- AddMetaData(seurat_obj, barcode_cluster, col.name = "coral_ground_truth_state")

  # --- 4. OPTIMIZED Heritable Gene Identification ---
  message("Step 4/6: Identifying heritable genes via Omega-squared...")
  cell_total <- ncol(sc_expression)
  mean_expression <- Matrix::rowMeans(sc_expression)
  SSt <- Matrix::rowSums((sc_expression - mean_expression)^2)
  lineage_number_actual <- length(unique(barcode_used))

  barcode_indices <- split(seq_along(barcode_used), barcode_used)
  lineage_means <- vapply(
    barcode_indices, function(idx) Matrix::rowMeans(sc_expression[, idx, drop = FALSE]),
    numeric(nrow(sc_expression))
  )
  lineage_mean_expression <- lineage_means[, match(as.character(barcode_used), names(barcode_indices))]
  SSw <- Matrix::rowSums((sc_expression - lineage_mean_expression)^2)

  SS <- data.frame(SSt = SSt, SSw = SSw, SSb = SSt - SSw, Mean = mean_expression, name = rownames(sc_expression))
  SS$Omega_square <- (SS$SSb / (cell_total - lineage_number_actual)) / (SS$SSt / (cell_total - 1))
  SS$CV <- sqrt(SS$SSt / SS$Mean) / sqrt(cell_total)
  SS_filter <- subset(SS, Mean > min_gene_mean)
  
  message("  - Running parallel permutation test for Omega-squared threshold...")
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  omega_square_fake_list <- foreach::foreach(i = 1:permutation_repeats, .packages="Matrix") %dopar% {
      barcode_fake <- sample(barcode_used)
      fake_indices <- split(seq_along(barcode_fake), barcode_fake)
      fake_means <- vapply(
        fake_indices, function(idx) Matrix::rowMeans(sc_expression[, idx, drop = FALSE]),
        numeric(nrow(sc_expression))
      )
      fake_mean_expr <- fake_means[, match(as.character(barcode_fake), names(fake_indices))]
      SSw_fake <- Matrix::rowSums((sc_expression - fake_mean_expr)^2)
      SSb_fake <- SSt - SSw_fake
      (SSb_fake / (cell_total - lineage_number_actual)) / (SSt / (cell_total - 1))
  }
  parallel::stopCluster(cl)
  omega_square_fake <- unlist(omega_square_fake_list)

  threshold <- 1.5 * stats::IQR(omega_square_fake, na.rm = TRUE) + stats::quantile(omega_square_fake, 0.75, na.rm = TRUE)

  # --- 5. Calculate MDS Coordinates ---
  message("Step 5/6: Calculating MDS coordinates...")
  mds_coords <- as.data.frame(stats::cmdscale(clone_dist_E, k = 10))
  colnames(mds_coords) <- paste0("MDS", 1:10)
  mds_coords$BarcodeID <- as.character(rownames(mds_coords))

  # --- 6. Store Results ---
  message("Step 6/6: Storing results in Seurat object misc slot...")
  ground_truth_results <- list(
    clone_energy_distance = clone_dist_E,
    clone_statistics = barcode_sta_df_filtered,
    mds_coordinates = mds_coords,
    heritable_genes_df = SS_filter,
    heritable_genes_threshold = threshold,
    omega_squared_null_distribution = omega_square_fake,
    parameters = list(
        clone_size_cutoff = clone_size_cutoff, num_states = num_states,
        hclust_method = hclust_method
    ),
    internal_data = list(
        sc_expression_subset = sc_expression,
        barcode_used = barcode_used
    )
  )
  seurat_obj@misc$CORAL_ground_truth_analysis <- ground_truth_results

  message("CORAL ground truth analysis complete. Results stored in 'seurat_obj@misc$CORAL_ground_truth_analysis'.")
  return(seurat_obj)
}

#' @title Analyze Gene Fluctuation Mode (Highly Optimized)
#' @description This optimized version performs hierarchical clustering only once,
#' then uses cutree() to efficiently get all cluster assignments for calculating
#' gene-specific fluctuation metrics (`omega_area`).
#'
#' @param seurat_obj A Seurat object that has been processed by `run_coral_ground_truth_analysis`.
#' @param n_cores An integer specifying the number of cores to use for parallel processing.
#'
#' @return An updated Seurat object with fluctuation analysis results.
#'
#' @export
analyze_gene_fluctuation <- function(seurat_obj, n_cores = NULL) {
    if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis)) {
        stop("Please run 'run_coral_ground_truth_analysis' first.")
    }
    message("Starting advanced analysis: Gene Fluctuation Mode (Optimized Version)...")

    res <- seurat_obj@misc$CORAL_ground_truth_analysis
    clone_dist_E <- res$clone_energy_distance
    lineage_total <- nrow(clone_dist_E)

    heritable_genes_df <- res$heritable_genes_df
    high_F_gene_names <- heritable_genes_df$name[heritable_genes_df$Omega_square > res$heritable_genes_threshold]

    if(length(high_F_gene_names) == 0) {
        warning("No heritable genes found above the threshold. Skipping fluctuation analysis.")
        return(seurat_obj)
    }

    sc_expression_ANOVA <- res$internal_data$sc_expression_subset[high_F_gene_names, ]
    barcode_used <- res$internal_data$barcode_used
    barcode_sta_df <- res$clone_statistics
    cell_total <- length(barcode_used)
    n_genes_anova <- nrow(sc_expression_ANOVA)

    message("  - Performing hierarchical clustering once...")
    hclust_result <- stats::hclust(as.dist(clone_dist_E), method = res$parameters$hclust_method)

    message("  - Pre-calculating all cluster partitions for gene analysis...")
    all_partitions_genes <- stats::cutree(hclust_result, k = 1:(lineage_total - 1))
    
    clone_id_map <- match(as.character(barcode_sta_df$BarcodeID), hclust_result$labels)
    all_partitions_reordered <- all_partitions_genes[clone_id_map, ]
    
    cell_to_clone_map <- match(barcode_used, barcode_sta_df$BarcodeID)
    
    if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    message("  - Calculating SSw across hierarchy levels (lightweight parallel loop)...")
    SSw_results_list <- foreach::foreach(
        j = 2:(lineage_total - 1),
        .packages = "Matrix"
    ) %dopar% {
        clone_clusters_j <- all_partitions_reordered[, j-1]
        barcode_cluster_j <- clone_clusters_j[cell_to_clone_map]

        valid_clusters <- sort(unique(barcode_cluster_j))
        cluster_means <- vapply(valid_clusters, function(k) {
            Matrix::rowMeans(sc_expression_ANOVA[, barcode_cluster_j == k, drop = FALSE])
        }, numeric(n_genes_anova))
        
        cluster_map <- match(barcode_cluster_j, valid_clusters)
        mean_expr_mapped <- cluster_means[, cluster_map]
        
        Matrix::rowSums((sc_expression_ANOVA - mean_expr_mapped)^2)
    }
    parallel::stopCluster(cl)

    SSw_matrix <- do.call(cbind, SSw_results_list)

    SSt_anova <- heritable_genes_df[high_F_gene_names, "SSt"]
    SSw_anova_last <- heritable_genes_df[high_F_gene_names, "SSw"]
    SSw_full_matrix <- cbind(SSt_anova, SSw_matrix, SSw_anova_last)

    lineage_numbers <- 1:lineage_total
    
    term1 <- SSt_anova %o% (cell_total - lineage_numbers)
    term2 <- cell_total * SSw_full_matrix
    gene_omega_square_matrix <- (term1 - term2) / (term1 + term2)

    norm_gene_omega_square <- gene_omega_square_matrix / gene_omega_square_matrix[, lineage_total]
    
    ref_curve <- (0:(lineage_total - 1)) / (lineage_total - 1)
    gene_omega_area <- rowSums(norm_gene_omega_square, na.rm = TRUE) - sum(ref_curve, na.rm = TRUE)
    gene_omega_area_norm <- gene_omega_area / lineage_total

    seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$omega_area <- NA
    match_indices <- match(high_F_gene_names, seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$name)
    seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$omega_area[match_indices] <- gene_omega_area_norm
    
    seurat_obj@misc$CORAL_ground_truth_analysis$fluctuation_analysis <- list(
        gene_omega_square_matrix = gene_omega_square_matrix,
        normalized_gene_omega_square = norm_gene_omega_square
    )

    message("Gene fluctuation analysis complete.")
    return(seurat_obj)
}

#' @title Visualize the Clone Energy Distance Heatmap
#' @description Creates a heatmap of the clone-clone energy distance matrix, visualizing
#'   the relationships between clones and the defined CORAL states.
#'
#' @param seurat_obj A Seurat object processed by `run_coral_ground_truth_analysis`.
#' @param ... Additional parameters passed to `ComplexHeatmap::Heatmap` for customization.
#'
#' @return A `Heatmap` object from the `ComplexHeatmap` package.
#'
#' @seealso \code{\link{run_coral_ground_truth_analysis}}, \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @export
visualize_clone_distance_heatmap <- function(seurat_obj, ...) {
    if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis)) {stop("Please run 'run_coral_ground_truth_analysis' first.")}
    res <- seurat_obj@misc$CORAL_ground_truth_analysis
    clone_dist_E <- res$clone_energy_distance
    num_states <- res$parameters$num_states
    col_fun <- circlize::colorRamp2(c(0, stats::median(clone_dist_E), max(clone_dist_E)), c("red", "white", "darkslateblue"))
    ht <- ComplexHeatmap::Heatmap(clone_dist_E, name = "E-distance", col = col_fun, row_split = num_states, column_split = num_states,
        show_row_names = FALSE, show_column_names = FALSE, row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        clustering_method_rows = res$parameters$hclust_method, clustering_method_columns = res$parameters$hclust_method, ...)
    return(ht)
}


#' @title Visualize the MDS Embedding of Clones
#' @description Generates a 2D scatter plot of clones based on the Multi-Dimensional Scaling (MDS)
#'   of the energy distance matrix.
#'
#' @param seurat_obj A Seurat object processed by `run_coral_ground_truth_analysis`.
#' @param color_by A string specifying how to color the points. Can be `"coral_state"` (default)
#'   or any valid column name from the Seurat object's metadata.
#' @param mds_dims A numeric vector of length 2 specifying which MDS dimensions to plot. Defaults to `c(1, 2)`.
#' @param pt_size The size of the points in the plot. Defaults to `2`.
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{run_coral_ground_truth_analysis}}
#'
#' @export
visualize_clone_mds <- function(seurat_obj, color_by = "coral_state", mds_dims = c(1, 2), pt_size = 2) {
    if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis)) { stop("Please run 'run_coral_ground_truth_analysis' first.") }
    res <- seurat_obj@misc$CORAL_ground_truth_analysis
    mds_df <- res$mds_coordinates
    color_label <- color_by
    plot_data <- mds_df
    if (color_by == "coral_state") {
        clone_to_state_map <- unique(seurat_obj@meta.data[, c("coral_barcode_numeric", "coral_ground_truth_state")])
        clone_to_state_map <- clone_to_state_map[clone_to_state_map$coral_ground_truth_state != 0, ]
        plot_data <- merge(plot_data, clone_to_state_map, by.x="BarcodeID", by.y="coral_barcode_numeric", all.x=TRUE)
        plot_data$Color <- as.factor(plot_data$coral_ground_truth_state)
        color_label <- "CORAL State"
        num_colors <- length(unique(stats::na.omit(plot_data$Color)))
        color_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_colors)
        p <- ggplot(plot_data, aes(x = .data[[paste0("MDS", mds_dims[1])]], y = .data[[paste0("MDS", mds_dims[2])]], color = .data$Color)) +
             scale_color_manual(values = color_palette, na.value="grey")
    } else if (color_by %in% colnames(seurat_obj@meta.data)) {
        clone_fate <- table(seurat_obj$coral_barcode_numeric, seurat_obj@meta.data[[color_by]])
        clone_fate_prop <- as.data.frame.matrix(clone_fate / rowSums(clone_fate))
        clone_fate_prop$BarcodeID <- rownames(clone_fate_prop)
        plot_data <- merge(plot_data, clone_fate_prop, by="BarcodeID", all.x=TRUE)
        first_level <- colnames(clone_fate_prop)[1]
        plot_data$Color <- plot_data[[first_level]]
        color_label <- paste("Proportion", first_level)
        p <- ggplot(plot_data, aes(x = .data[[paste0("MDS", mds_dims[1])]], y = .data[[paste0("MDS", mds_dims[2])]], color = .data$Color)) +
             scale_color_distiller(palette = "RdYlBu")
    } else { stop("'color_by' must be 'coral_state' or a valid column name.") }
    p <- p + geom_point(size = pt_size) + labs(x = paste("MDS", mds_dims[1]), y = paste("MDS", mds_dims[2]), color = color_label,
       title = "MDS Embedding of Clones", subtitle = paste("Colored by", color_by)) + theme_classic()
    return(p)
}


#' @title Plot Heritable Gene Effect Size Distribution
#' @description Creates a density plot comparing the distribution of observed Omega-squared values
#'   for all genes against the null distribution generated from permutations.
#'
#' @param seurat_obj A Seurat object processed by `run_coral_ground_truth_analysis`.
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{run_coral_ground_truth_analysis}}
#'
#' @export
plot_heritable_gene_distribution <- function(seurat_obj) {
    if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis)) { stop("Please run 'run_coral_ground_truth_analysis' first.") }
    results <- seurat_obj@misc$CORAL_ground_truth_analysis
    scores_df <- results$heritable_genes_df
    null_dist <- results$omega_squared_null_distribution
    threshold <- results$heritable_genes_threshold
    plot_df <- data.frame(omega_square = c(scores_df$Omega_square, null_dist),
        source = c(rep("Observed", nrow(scores_df)), rep("Permuted Null", length(null_dist))))
    plot_df <- plot_df[!is.na(plot_df$omega_square), ]
    p <- ggplot(plot_df, aes(x = .data$omega_square, fill = .data$source, color = .data$source)) +
        geom_density(alpha = 0.5, adjust = 2) +
        geom_vline(xintercept = threshold, linetype = "dashed", color = "red", size = 1) +
        scale_x_log10(limits = c(NA, 1), name = "Omega-squared") +
        scale_fill_manual(values = c("Observed" = "coral", "Permuted Null" = "grey60")) +
        scale_color_manual(values = c("Observed" = "coral", "Permuted Null" = "grey60")) +
        labs(title = "Heritable Gene Effect Size Distribution", subtitle = "Red dashed line indicates significance threshold", y = "Gene Density") +
        theme_classic() + theme(legend.title = element_blank())
    return(p)
}


#' @title Visualize CORAL States on UMAP
#' @description Overlays the calculated CORAL states on a UMAP projection of the cells.
#'   Note: Requires UMAP to be computed on the Seurat object first.
#'
#' @param seurat_obj A Seurat object processed by `run_coral_ground_truth_analysis`.
#' @param ... Additional parameters passed to `Seurat::DimPlot`.
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link[Seurat]{DimPlot}}, \code{\link{run_coral_ground_truth_analysis}}
#'
#' @export
visualize_coral_states_umap <- function(seurat_obj, ...) {
    if (!"coral_ground_truth_state" %in% colnames(seurat_obj@meta.data)) { stop("Please run 'run_coral_ground_truth_analysis' first.") }
    if (!"umap" %in% names(seurat_obj@reductions)) { warning("UMAP reduction not found. Please run RunUMAP() on the Seurat object first.")}
    Seurat::DimPlot(seurat_obj, group.by = "coral_ground_truth_state", label = TRUE, ...) +
        ggplot2::ggtitle("CORAL States on UMAP")
}


#' @title Plot Confusion Matrix of CORAL States vs. Cell Types
#' @description Generates a heatmap showing the overlap (confusion matrix) between the
#'   computed CORAL states and another cell annotation, such as cell type.
#'
#' @param seurat_obj A Seurat object processed by `run_coral_ground_truth_analysis`.
#' @param celltype_col A string specifying the metadata column with the cell type or other annotations to compare against.
#' @param ... Additional parameters passed to `ComplexHeatmap::pheatmap`.
#'
#' @return A `Heatmap` object from the `ComplexHeatmap` package.
#'
#' @seealso \code{\link[ComplexHeatmap]{pheatmap}}
#'
#' @export
plot_state_celltype_confusion <- function(seurat_obj, celltype_col, ...) {
    if (!"coral_ground_truth_state" %in% colnames(seurat_obj@meta.data) || !celltype_col %in% colnames(seurat_obj@meta.data)) {
        stop("Ensure both 'coral_ground_truth_state' and 'celltype_col' exist in metadata.")}
    plot_df <- seurat_obj@meta.data[seurat_obj$coral_ground_truth_state != 0, ]
    confusion_matrix <- as.matrix(table(plot_df$coral_ground_truth_state, plot_df[[celltype_col]]))
    ComplexHeatmap::pheatmap(confusion_matrix, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, ...)
}


#' @title Analyze Gene Fluctuation Mode (Highly Optimized)
#' @description This optimized version performs hierarchical clustering only once,
#' then uses cutree() to efficiently get all cluster assignments for calculating
#' gene-specific fluctuation metrics (`omega_area`).
#'
#' @param seurat_obj A Seurat object that has been processed by `run_coral_ground_truth_analysis`.
#' @param n_cores An integer specifying the number of cores to use for parallel processing.
#'
#' @return An updated Seurat object with fluctuation analysis results.
#'
#' @export
analyze_gene_fluctuation <- function(seurat_obj, n_cores = NULL) {
    if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis)) {
        stop("Please run 'run_coral_ground_truth_analysis' first.")
    }
    message("Starting advanced analysis: Gene Fluctuation Mode (Optimized Version)...")

    res <- seurat_obj@misc$CORAL_ground_truth_analysis
    clone_dist_E <- res$clone_energy_distance
    lineage_total <- nrow(clone_dist_E)

    heritable_genes_df <- res$heritable_genes_df
    high_F_gene_names <- heritable_genes_df$name[heritable_genes_df$Omega_square > res$heritable_genes_threshold]

    if(length(high_F_gene_names) == 0) {
        warning("No heritable genes found above the threshold. Skipping fluctuation analysis.")
        return(seurat_obj)
    }

    sc_expression_ANOVA <- res$internal_data$sc_expression_subset[high_F_gene_names, ]
    barcode_used <- res$internal_data$barcode_used
    barcode_sta_df <- res$clone_statistics
    cell_total <- length(barcode_used)
    n_genes_anova <- nrow(sc_expression_ANOVA)

    message("  - Performing hierarchical clustering once...")
    hclust_result <- stats::hclust(as.dist(clone_dist_E), method = res$parameters$hclust_method)

    message("  - Pre-calculating all cluster partitions for gene analysis...")
    all_partitions_genes <- stats::cutree(hclust_result, k = 1:(lineage_total - 1))
    
    clone_id_map <- match(as.character(barcode_sta_df$BarcodeID), hclust_result$labels)
    all_partitions_reordered <- all_partitions_genes[clone_id_map, ]
    
    cell_to_clone_map <- match(barcode_used, barcode_sta_df$BarcodeID)
    
    if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    message("  - Calculating SSw across hierarchy levels (lightweight parallel loop)...")
    SSw_results_list <- foreach::foreach(
        j = 2:(lineage_total - 1),
        .packages = "Matrix"
    ) %dopar% {
        clone_clusters_j <- all_partitions_reordered[, j-1]
        barcode_cluster_j <- clone_clusters_j[cell_to_clone_map]

        valid_clusters <- sort(unique(barcode_cluster_j))
        cluster_means <- vapply(valid_clusters, function(k) {
            Matrix::rowMeans(sc_expression_ANOVA[, barcode_cluster_j == k, drop = FALSE])
        }, numeric(n_genes_anova))
        
        cluster_map <- match(barcode_cluster_j, valid_clusters)
        mean_expr_mapped <- cluster_means[, cluster_map]
        
        Matrix::rowSums((sc_expression_ANOVA - mean_expr_mapped)^2)
    }
    parallel::stopCluster(cl)

    SSw_matrix <- do.call(cbind, SSw_results_list)

    SSt_anova <- heritable_genes_df[high_F_gene_names, "SSt"]
    SSw_anova_last <- heritable_genes_df[high_F_gene_names, "SSw"]
    SSw_full_matrix <- cbind(SSt_anova, SSw_matrix, SSw_anova_last)

    lineage_numbers <- 1:lineage_total
    
    term1 <- SSt_anova %o% (cell_total - lineage_numbers)
    term2 <- cell_total * SSw_full_matrix
    gene_omega_square_matrix <- (term1 - term2) / (term1 + term2)

    norm_gene_omega_square <- gene_omega_square_matrix / gene_omega_square_matrix[, lineage_total]
    
    ref_curve <- (0:(lineage_total - 1)) / (lineage_total - 1)
    gene_omega_area <- rowSums(norm_gene_omega_square, na.rm = TRUE) - sum(ref_curve, na.rm = TRUE)
    gene_omega_area_norm <- gene_omega_area / lineage_total

    seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$omega_area <- NA
    match_indices <- match(high_F_gene_names, seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$name)
    seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$omega_area[match_indices] <- gene_omega_area_norm
    
    seurat_obj@misc$CORAL_ground_truth_analysis$fluctuation_analysis <- list(
        gene_omega_square_matrix = gene_omega_square_matrix,
        normalized_gene_omega_square = norm_gene_omega_square
    )

    message("Gene fluctuation analysis complete.")
    return(seurat_obj)
}

#' @title Visualize Gene Expression on Clone MDS Plot
#' @description Overlays the pseudobulk expression of a specific gene onto the
#' MDS embedding of clones. The required pseudobulk matrix is now calculated internally.
#'
#' @param seurat_obj A Seurat object processed by the CORAL analysis.
#' @param gene A string specifying the name of the gene to visualize.
#' @param log_transform Logical, whether to log2(x+1) transform the expression
#' values for visualization. Defaults to TRUE.
#' @param mds_dims A numeric vector of length 2 indicating which MDS dimensions to plot.
#' Defaults to `c(1, 2)`.
#' @param pt_size The size of the points in the plot.
#'
#' @return A ggplot object.
#' @export
visualize_gene_mds <- function(seurat_obj,
                               gene,
                               log_transform = TRUE,
                               mds_dims = c(1, 2),
                               pt_size = 2) {

  if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis)) {
    stop("CORAL analysis results not found. Please run run_coral_ground_truth_analysis() first.")
  }

  # --- 内部集成：计算并修正伪批量矩阵 ---
  results <- seurat_obj@misc$CORAL_ground_truth_analysis
  sc_expression_subset <- results$internal_data$sc_expression_subset
  barcode_used <- results$internal_data$barcode_used
  
  temp_seurat_for_agg <- CreateSeuratObject(counts = sc_expression_subset)
  temp_seurat_for_agg <- AddMetaData(temp_seurat_for_agg, metadata = as.character(barcode_used), col.name = "clone_id")
  Idents(temp_seurat_for_agg) <- "clone_id"
  
  pseudobulk_matrix <- AggregateExpression(temp_seurat_for_agg, assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA
  colnames(pseudobulk_matrix) <- gsub("^g", "", colnames(pseudobulk_matrix))
  # --- 逻辑结束 ---

  if (!gene %in% rownames(pseudobulk_matrix)) {
    stop(paste("Gene '", gene, "' not found in the expression data of the analyzed clones."))
  }

  mds_df <- results$mds_coordinates
  
  if (!"BarcodeID" %in% colnames(mds_df)) {
    stop("'BarcodeID' column not found in MDS coordinates data frame.")
  }

  common_clones <- intersect(mds_df$BarcodeID, colnames(pseudobulk_matrix))
  if (length(common_clones) == 0) {
    stop("No common clone IDs found between MDS coordinates and pseudobulk matrix after processing.")
  }

  mds_df_filtered <- mds_df[mds_df$BarcodeID %in% common_clones, ]
  expression_values <- pseudobulk_matrix[gene, mds_df_filtered$BarcodeID]

  if (log_transform) {
    expression_values <- log2(expression_values + 1)
  }

  mds_df_filtered$expression <- expression_values
  
  dim1_name <- paste0("MDS", mds_dims[1])
  dim2_name <- paste0("MDS", mds_dims[2])
  colnames(mds_df_filtered)[1:2] <- c(dim1_name, dim2_name)

  p <- ggplot(mds_df_filtered, aes(x = .data[[dim1_name]], y = .data[[dim2_name]])) +
    geom_point(aes(color = expression), size = pt_size) +
    scale_colour_distiller(palette = "Spectral") +
    labs(
      title = paste("MDS Embedding Colored by", gene, "Expression"),
      x = paste("MDS", mds_dims[1]),
      y = paste("MDS", mds_dims[2]),
      color = ifelse(log_transform, "log2(Expr+1)", "Expression")
    ) +
    theme_classic()

  return(p)
}

#' @title Plot Single Gene Omega-squared Flucuation Curve
#' @description Plots the normalized Omega-squared curve for a single heritable gene
#' across different levels of clone clustering, comparing it to a reference line.
#'
#' @param seurat_obj A Seurat object processed by the CORAL analysis, including
#' the advanced fluctuation analysis.
#' @param gene A string specifying the name of the heritable gene to visualize.
#'
#' @return A ggplot object.
#'
#' @export
plot_gene_omega_curve <- function(seurat_obj, gene) {

  # --- FIX 1: Corrected data path and object name ---
  # Check for the fluctuation_analysis list and the correct matrix within it.
  if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis$fluctuation_analysis$normalized_gene_omega_square)) {
    stop("Fluctuation analysis results not found. Please run analyze_gene_fluctuation() first.")
  }

  norm_omega_matrix <- seurat_obj@misc$CORAL_ground_truth_analysis$fluctuation_analysis$normalized_gene_omega_square

  # --- FIX 2: Check against rownames, not colnames ---
  if (!gene %in% rownames(norm_omega_matrix)) {
    stop(paste("Gene '", gene, "' not found in the fluctuation analysis results.",
               "This plot is only for genes identified as heritable."))
  }

  # --- FIX 3: Rebuild the plotting data frame from the matrix structure ---
  lineage_total <- ncol(norm_omega_matrix)
  
  plot_df <- data.frame(
    # The x-axis is the number of clusters, from 1 to the total number of clones
    cluster_num = 1:lineage_total,
    # The gene's value is the corresponding row from the matrix
    gene_value = norm_omega_matrix[gene, ],
    # The reference line is a simple diagonal from 0 to 1
    ref_value = (0:(lineage_total - 1)) / (lineage_total - 1)
  )

  # The ggplot call is now based on the correctly structured data frame
  p <- ggplot(plot_df, aes(x = .data$cluster_num)) +
    geom_line(aes(y = .data$gene_value, color = "Gene"), size = 1.2) +
    geom_line(aes(y = .data$ref_value, color = "Reference"), linetype = "dashed", size = 1.2) +
    scale_color_manual(
      name = "Trace",
      values = c("Gene" = "cyan4", "Reference" = "darkgray")
    ) +
    labs(
      title = paste("Omega-squared Fluctuation for", gene),
      x = "Number of Lineage Clusters",
      y = "Normalized Omega-squared"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 15),
      axis.line = element_line(size = 1),
      legend.position = "top"
    )

  return(p)
}

#' @title Create a Comprehensive Gene Analysis Dashboard
#' @description Combines multiple plots for a single gene into one comprehensive view.
#' It now includes an AUC value annotation on the fluctuation curve plot.
#' The required pseudobulk matrix is calculated internally.
#'
#' @param seurat_obj A Seurat object processed by the CORAL analysis.
#' @param gene A string specifying the name of the gene to analyze.
#' @param num_ridge_idents An integer for the number of top clones to show in the
#' ridge plot. Defaults to 12.
#'
#' @return A plot object created by `ggpubr::ggarrange`.
#' @export
plot_gene_dashboard <- function(seurat_obj,
                                gene,
                                num_ridge_idents = 12) {

  results <- seurat_obj@misc$CORAL_ground_truth_analysis
  
  # --- Internal logic: Calculate and clean the pseudobulk matrix ---
  sc_expression_subset <- results$internal_data$sc_expression_subset
  barcode_used <- results$internal_data$barcode_used
  
  temp_seurat_for_agg <- CreateSeuratObject(counts = sc_expression_subset)
  temp_seurat_for_agg <- AddMetaData(temp_seurat_for_agg, metadata = as.character(barcode_used), col.name = "clone_id")
  Idents(temp_seurat_for_agg) <- "clone_id"
  
  pseudobulk_matrix <- AggregateExpression(temp_seurat_for_agg, assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA
  colnames(pseudobulk_matrix) <- gsub("^g", "", colnames(pseudobulk_matrix))
  # --- End of internal logic ---

  # 1. Generate Omega Curve Plot
  p_omega <- plot_gene_omega_curve(seurat_obj, gene) +
    theme(legend.position = "none")

  # --- [FEATURE]: Extract AUC value and add annotation to p_omega plot ---
  heritable_genes_df <- results$heritable_genes_df
  # Ensure the corresponding gene info exists in heritable_genes_df
  if (gene %in% heritable_genes_df$name) {
    gene_info <- heritable_genes_df[heritable_genes_df$name == gene, ]
    auc_value <- gene_info$omega_area
    
    # Check if the AUC value exists and is not NA
    if (length(auc_value) == 1 && !is.na(auc_value)) {
      # Create the label text
      auc_text <- sprintf("AUC: %.3f", auc_value)
      
      # --- [MODIFIED]: Use annotate() to add the text to the BOTTOM-RIGHT corner of the plot ---
      p_omega <- p_omega +
        ggplot2::annotate("text", x = Inf, y = -Inf, label = auc_text, # MODIFIED: y = -Inf to anchor at the bottom
                          hjust = 1.1, vjust = -0.5, size = 4.5, color = "black", # MODIFIED: vjust to push text UP from the bottom
                          fontface = "bold")
    }
  }
  # --- End of feature section ---

  # 2. Generate UMAP Feature Plot
  p_umap <- FeaturePlot(seurat_obj, features = gene, pt.size = 1, order = TRUE) +
    NoAxes() +
    ggtitle(NULL)

  # 3. Generate MDS Expression Plot
  p_mds <- visualize_gene_mds(seurat_obj, gene = gene) +
    ggtitle(NULL)

  # 4. Generate Ridge Plot for top clones
  clone_stats <- results$clone_statistics
  top_clones <- head(clone_stats$BarcodeID, num_ridge_idents)
  
  Idents(seurat_obj) <- "coral_barcode_numeric"
  p_ridge <- RidgePlot(seurat_obj, features = gene, idents = as.character(top_clones),
                       sort = "decreasing") +
    NoLegend() +
    ggtitle("Expression in Top Clones") +
    theme(axis.title.y = element_blank())

  # Arrange all plots into a dashboard (panel labels removed)
  dashboard <- ggpubr::ggarrange(
    p_omega, p_umap, p_mds, p_ridge,
    ncol = 2, nrow = 2
  )

  # Add a main title
  final_plot <- ggpubr::annotate_figure(
    dashboard,
    top = ggpubr::text_grob(paste("Comprehensive Analysis of Gene:", gene),
                           face = "bold", size = 16)
  )

  return(final_plot)
}

#' @title Plot Gene Fluctuation Mode
#' @description Creates a scatter plot visualizing the gene fluctuation mode.
#' This version includes a reference line for the significance threshold (y-axis).
#'
#' @param seurat_obj A Seurat object that has been run through `analyze_gene_fluctuation`.
#' @param n_genes_to_show An integer. Number of top genes to display.
#' @param genes_to_highlight A character vector of specific gene names to label.
#'
#' @return A ggplot object.
#' @export
plot_gene_fluctuation_mode <- function(seurat_obj, n_genes_to_show = 1000, genes_to_highlight = NULL) {
    if (is.null(seurat_obj@misc$CORAL_ground_truth_analysis$heritable_genes_df$omega_area)) {
        stop("Please run 'analyze_gene_fluctuation' first.")
    }
    
    results <- seurat_obj@misc$CORAL_ground_truth_analysis
    plot_df <- results$heritable_genes_df
    plot_df <- plot_df[!is.na(plot_df$omega_area), ]
    plot_df <- plot_df[order(-plot_df$Omega_square), ]
    
    n_genes_to_show <- min(n_genes_to_show, nrow(plot_df))
    plot_df_subset <- plot_df[1:n_genes_to_show, ]
    
    label_df <- plot_df[plot_df$name %in% genes_to_highlight, ]

    threshold <- results$heritable_genes_threshold

    p <- ggplot(plot_df_subset, aes(x = .data$omega_area, y = .data$Omega_square)) +
        # Only the horizontal line for the heritability threshold remains
        geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 0.8) +
        
        # The geom_vline for the global reference has been removed
        
        # Scatter plot layers
        geom_point(color = "grey30", size = 2, alpha = 0.8) +
        geom_point(data = label_df, color = "salmon", size = 3) +
        ggrepel::geom_label_repel(data = label_df, aes(label = .data$name),
                                  max.overlaps = Inf, color = "salmon4", size = 4) +
        theme_classic() +
        labs(x = "Normalized AUC Shape (omega_area)",
             y = "Lineage Effect Size (Omega-squared)",
             title = "Gene Fluctuation Mode",
             subtitle = "Red dashed line indicates heritability threshold") # Updated subtitle
    
    return(p)
}