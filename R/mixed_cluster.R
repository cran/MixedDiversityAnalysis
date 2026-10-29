
#' Perform Hierarchical Clustering on Mixed-Trait Data
#'
#' Performs hierarchical clustering using a mixed-distance object
#' created by mixed_distance().
#'
#' @param distance_object An object of class "mixed_distance".
#' @param k Number of clusters, or "optimal" to select the number
#'   using average silhouette width.
#' @param method Hierarchical linkage method: "average", "complete",
#'   "single", or "mcquitty".
#' @param k_range Integer vector of candidate cluster numbers when
#'   k = "optimal". If NULL, a suitable range is generated.
#' @param min_cluster_size Minimum permitted cluster size during
#'   automatic selection. Defaults to 1.
#'
#' @return An object of class "mixed_cluster" containing the
#'   clustering tree, cluster assignments, cluster sizes, selected
#'   number of clusters, linkage method, silhouette results,
#'   comparison of candidate solutions, and minimum cluster size.
#'
#' @importFrom stats hclust cutree
#' @importFrom cluster silhouette
#'
#' @examples
#' # Prepare the same rice dataset used in mixed_distance()
#' rice_data <- data.frame(
#'   PlantHeight = c(
#'     90, 100, 110, 95, 85, 105, 115, 92, 98, 108,
#'     88, 102, 112, 97, 83, 107, 118, 93, 101, 99
#'   ),
#'   GrainYield = c(
#'     25, 30, 35, 28, 22, 32, 38, 26, 29, 34,
#'     24, 31, 36, 27, 21, 33, 40, 25, 30, 28
#'   ),
#'   GrainType = factor(c(
#'     "Long", "Long", "Bold", "Bold", "Medium",
#'     "Long", "Bold", "Medium", "Long", "Bold",
#'     "Medium", "Long", "Bold", "Medium", "Bold",
#'     "Long", "Bold", "Medium", "Long", "Bold"
#'   )),
#'   DiseaseResponse = factor(c(
#'     "Resistant", "Susceptible", "Resistant", "Susceptible",
#'     "Resistant", "Resistant", "Susceptible", "Resistant",
#'     "Susceptible", "Resistant", "Susceptible", "Resistant",
#'     "Resistant", "Susceptible", "Resistant", "Susceptible",
#'     "Resistant", "Susceptible", "Resistant", "Susceptible"
#'   )),
#'   Maturity = ordered(
#'     c(
#'       "Early", "Medium", "Late", "Medium", "Early",
#'       "Late", "Late", "Early", "Medium", "Late",
#'       "Early", "Medium", "Late", "Medium", "Early",
#'       "Late", "Late", "Early", "Medium", "Medium"
#'     ),
#'     levels = c("Early", "Medium", "Late")
#'   )
#' )
#'
#' rownames(rice_data) <- paste0("Rice_", sprintf("%02d", 1:20))
#'
#' trait_types <- c(
#'   "quantitative", "quantitative", "qualitative",
#'   "binary", "ordinal"
#' )
#'
#' # Step 1: Calculate dissimilarities
#' dist_result <- mixed_distance(rice_data, trait_types)
#'
#' # Step 2: Cluster the rice genotypes
#' cluster_result <- mixed_cluster(
#'   dist_result,
#'   k = 3,
#'   method = "average"
#' )
#'
#' # Step 3: Inspect cluster assignments
#' print(cluster_result$clusters)


#' @export
mixed_cluster <- function(
    distance_object,
    k = "optimal",
    method = "average",
    k_range = NULL,
    min_cluster_size = 1
) {
  
  # 1. Validate distance object
  if (!inherits(distance_object, "mixed_distance")) {
    stop("'distance_object' must be created by mixed_distance().")
  }
  
  if (is.null(distance_object$dissimilarity)) {
    stop("The distance object has no dissimilarity matrix.")
  }
  
  dmat <- distance_object$dissimilarity
  
  if (!is.matrix(dmat) || nrow(dmat) != ncol(dmat)) {
    stop("The dissimilarity matrix must be square.")
  }
  
  n <- nrow(dmat)
  
  if (n < 3) {
    stop("At least three observations are required for clustering.")
  }
  
  if (anyNA(dmat) || any(!is.finite(dmat))) {
    stop(
      "The dissimilarity matrix contains missing or non-finite values."
    )
  }
  
  if (any(dmat < 0)) {
    stop("Dissimilarities must be non-negative.")
  }
  
  if (!isTRUE(all.equal(dmat, t(dmat)))) {
    stop("The dissimilarity matrix must be symmetric.")
  }
  
  
  if (any(abs(diag(dmat)) > 1e-10)) {
    stop("The diagonal of the dissimilarity matrix must be zero.")
  }
  
  
  if (any(dmat > 1)) {
    stop("Dissimilarities must not exceed 1.")
  }
  
  # 2. Validate linkage method
  allowed_methods <- c(
    "average",
    "complete",
    "single",
    "mcquitty"
  )
  
  if (length(method) != 1 ||
      is.na(method) ||
      !method %in% allowed_methods) {
    stop(
      "'method' must be average, complete, single, or mcquitty."
    )
  }
  
  # 3. Validate minimum cluster size
  if (length(min_cluster_size) != 1 ||
      !is.numeric(min_cluster_size) ||
      is.na(min_cluster_size) ||
      !is.finite(min_cluster_size) ||
      min_cluster_size < 1 ||
      min_cluster_size != floor(min_cluster_size)) {
    stop("'min_cluster_size' must be a positive integer.")
  }
  
  # 4. Convert matrix to a dist object for clustering
  distance <- stats::as.dist(dmat)
  
  tree <- stats::hclust(
    distance,
    method = method
  )
  
  # Helper to calculate silhouette summary
  evaluate_solution <- function(cluster_ids) {
    
    sizes <- table(cluster_ids)
    
    if (length(sizes) < 2 || length(sizes) >= n) {
      return(NULL)
    }
    
    if (min(sizes) < min_cluster_size) {
      return(NULL)
    }
    
    sil <- cluster::silhouette(
      cluster_ids,
      distance
    )
    
    list(
      silhouette = sil,
      silhouette_mean = mean(sil[, "sil_width"]),
      sizes = sizes
    )
  }
  
  # 5. Choose a fixed or optimal number of clusters
  if (identical(k, "optimal")) {
    
    if (is.null(k_range)) {
      k_range <- seq.int(2, n - 1)
    }
    
    if (!is.numeric(k_range) ||
        length(k_range) == 0 ||
        anyNA(k_range) ||
        any(!is.finite(k_range)) ||
        any(k_range != floor(k_range)) ||
        any(k_range < 2 | k_range >= n)) {
      stop("'k_range' must contain integers from 2 to n - 1.")
    }
    
    k_range <- sort(unique(as.integer(k_range)))
    
    comparisons <- data.frame(
      k = integer(0),
      silhouette_mean = numeric(0)
    )
    
    solutions <- list()
    
    for (candidate_k in k_range) {
      
      candidate_clusters <- stats::cutree(
        tree,
        k = candidate_k
      )
      
      evaluation <- evaluate_solution(candidate_clusters)
      
      if (!is.null(evaluation)) {
        comparisons <- rbind(
          comparisons,
          data.frame(
            k = candidate_k,
            silhouette_mean = evaluation$silhouette_mean
          )
        )
        
        solutions[[as.character(candidate_k)]] <- list(
          clusters = candidate_clusters,
          evaluation = evaluation
        )
      }
    }
    
    if (nrow(comparisons) == 0) {
      stop(
        "No candidate solution satisfies min_cluster_size. ",
        "Try a smaller minimum size or a different k_range."
      )
    }
    
    # In ties, choose the smaller number of clusters
    best_row <- order(
      -comparisons$silhouette_mean,
      comparisons$k
    )[1]
    
    selected_k <- comparisons$k[best_row]
    best_solution <- solutions[[as.character(selected_k)]]
    
    clusters <- best_solution$clusters
    silhouette_result <- best_solution$evaluation$silhouette
    silhouette_mean <- best_solution$evaluation$silhouette_mean
    
  } else {
    
    if (!is.numeric(k) ||
        length(k) != 1 ||
        is.na(k) ||
        !is.finite(k) ||
        k != floor(k) ||
        k < 2 ||
        k >= n) {
      stop("'k' must be an integer between 2 and n - 1, or 'optimal'.")
    }
    
    selected_k <- as.integer(k)
    
    clusters <- stats::cutree(
      tree,
      k = selected_k
    )
    
    evaluation <- evaluate_solution(clusters)
    
    if (is.null(evaluation)) {
      stop(
        "The selected clusters cannot be evaluated. ",
        "Check min_cluster_size and the number of clusters."
      )
    }
    
    silhouette_result <- evaluation$silhouette
    silhouette_mean <- evaluation$silhouette_mean
    
    comparisons <- data.frame(
      k = selected_k,
      silhouette_mean = silhouette_mean
    )
  }
  
  # 6. Calculate cluster sizes
  cluster_sizes <- table(clusters)
  
  # 7. Return structured result
  result <- list(
    tree = tree,
    clusters = clusters,
    sizes = cluster_sizes,
    k = selected_k,
    method = method,
    silhouette = silhouette_result,
    silhouette_mean = silhouette_mean,
    comparison = comparisons,
    min_cluster_size = min_cluster_size
  )
  
  class(result) <- "mixed_cluster"
  
  return(result)
}
