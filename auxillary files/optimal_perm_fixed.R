optimal_perm <- function(node.clusters.true, hyper.clusters.true, c_t, d_t, constants = constants2) {
  library(clue)  # For solve_LSAP()
  
  # Extract cluster counts
  K <- as.numeric(constants[4])  # number of node clusters
  G <- as.numeric(constants[5])  # number of hyperedge clusters
  
  # Normalize soft assignments row-wise
  e$c_t <- e$c_t / rowSums(e$c_t)
  e$d_t <- e$d_t / rowSums(e$d_t)

  # Get hard labels (most probable cluster)
  get_hard_labels <- function(prob_matrix) {
    numeric_mat <- matrix(as.numeric(prob_matrix), nrow = nrow(prob_matrix), ncol = ncol(prob_matrix))
    max.col(numeric_mat, ties.method = "first")
  }
  
  # Optimize label permutation using Hungarian algorithm
  optimize_labels <- function(true_labels, prob_matrix, num_clusters) {
    predicted_labels <- get_hard_labels(prob_matrix)
    
    # Ensure labels are not NA and matched in size
    valid_idx <- which(!is.na(true_labels) & !is.na(predicted_labels))
    true <- true_labels[valid_idx]
    pred <- predicted_labels[valid_idx]
    
    # Create confusion matrix
    confusion_mat <- table(true, pred)
    
    # Fill into a square matrix if needed (for LSAP to work)
    confusion_square <- matrix(0, nrow = num_clusters, ncol = num_clusters)
    rownames(confusion_square) <- 1:num_clusters
    colnames(confusion_square) <- 1:num_clusters
    confusion_square[rownames(confusion_mat), colnames(confusion_mat)] <- confusion_mat
    
    # Solve assignment problem: find best permutation
    matching <- solve_LSAP(confusion_square, maximum = TRUE)
    best_permutation <- as.integer(matching)
    trace_val <- sum(confusion_square[cbind(1:num_clusters, matching)])
    
    return(list(
      optimal_labels = best_permutation,
      max_trace = trace_val,
      confusion_mat = confusion_square
    ))
  }
  
  # Apply optimization to node and hyperedge clustering
  node_result <- optimize_labels(node.clusters.true, e$c_t, K)
  cat("Optimal node labels:", node_result$optimal_labels, "\n")
  cat("Max trace for nodes:", node_result$max_trace, "\n")
  
  hyper_result <- optimize_labels(hyper.clusters.true, e$d_t, G)
  cat("Optimal hyperedge labels:", hyper_result$optimal_labels, "\n")
  cat("Max trace for hyperedges:", hyper_result$max_trace, "\n")
  
  return(list(
    node_labels = node_result$optimal_labels,
    hyper_labels = hyper_result$optimal_labels,
    node_trace = node_result$max_trace,
    hyper_trace = hyper_result$max_trace,
    confusion_mat_hyper = hyper_result$confusion_mat,
    confusion_mat_node = node_result$confusion_mat
  ))
}
