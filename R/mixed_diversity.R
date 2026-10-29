#' Measure Trait Contributions to Cluster Separation
#'
#' Quantifies how strongly individual traits distinguish clusters
#' using separation scores appropriate for different trait types.
#'
#' @param data A data frame containing the traits to be evaluated.
#' @param types A character vector specifying the type of each trait.
#'   Supported types are "quantitative", "qualitative", "binary",
#'   and "ordinal".
#' @param clusters A vector assigning each observation to a cluster.
#'   Its length must equal the number of rows in `data`, and values
#'   must not be missing.
#'
#' @return A data frame with three columns:
#'   \describe{
#'     \item{trait}{The name of each trait.}
#'     \item{type}{The specified trait type.}
#'     \item{separation}{A descriptive separation score between zero
#'       and one, or `NA` when a score cannot be calculated for an
#'       entirely missing trait.}
#'   }
#'   Traits are ordered by decreasing separation score.
#'
#' @details
#' For quantitative and ordinal traits, separation is measured as the
#' proportion of total variation attributable to differences between
#' clusters, using the ratio of between-cluster sum of squares to
#' total sum of squares.
#'
#' For qualitative and binary traits, separation is measured using
#' squared Cramer's V, calculated from the Pearson chi-squared
#' statistic. This score describes the association between trait
#' categories and cluster membership.
#'
#' Scores range from zero to one, with larger values indicating
#' stronger observed separation among clusters. These scores are
#' descriptive measures of association and do not establish causality
#' or statistical significance.

#'
#' @examples
#' # Prepare a dataset of 20 rice genotypes
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
#' # Define the type of each trait
#' trait_types <- c(
#'   "quantitative", "quantitative", "qualitative",
#'   "binary", "ordinal"
#' )
#'
#' # Calculate mixed-trait dissimilarities
#' dist_result <- mixed_distance(rice_data, trait_types)
#'
#' # Identify three clusters
#' cluster_result <- mixed_cluster(
#'   dist_result,
#'   k = 3,
#'   method = "average"
#' )
#'
#' # Measure the separation contributed by each trait
#' trait_result <- mixed_diversity(
#'   data = rice_data,
#'   types = trait_types,
#'   clusters = cluster_result$clusters
#' )
#'
#' # Display trait separation scores
#' print(trait_result)

#' @export
mixed_diversity <- function(data, types, clusters) {
  
  # 1. Validate input data
  if (!is.data.frame(data) || nrow(data) < 3 || ncol(data) < 1) {
    stop("'data' must be a data frame with at least 3 rows and 1 column.")
  }
  
  if (!is.character(types) ||
      length(types) != ncol(data) ||
      anyNA(types)) {
    stop("'types' must specify the type of every column in data.")
  }
  
  allowed_types <- c(
    "quantitative", "qualitative", "binary", "ordinal"
  )
  
  if (any(!types %in% allowed_types)) {
    stop("Unsupported trait type in 'types'.")
  }
  
  if (length(clusters) != nrow(data) ||
      anyNA(clusters)) {
    stop("'clusters' must have one non-missing value per observation.")
  }
  
  clusters <- as.factor(clusters)
  
  if (nlevels(droplevels(clusters)) < 2) {
    stop("'clusters' must contain at least two clusters.")
  }
  
  if (any(table(clusters) < 1)) {
    stop("Each cluster must contain at least one observation.")
  }
  
  # 2. Calculate one separation score per trait
  scores <- numeric(ncol(data))
  
  for (j in seq_len(ncol(data))) {
    
    x <- data[[j]]
    trait_type <- types[j]
    
    if (all(is.na(x))) {
      scores[j] <- NA_real_
      next
    }
    
    if (trait_type == "quantitative") {
      
      if (!is.numeric(x)) {
        stop("Quantitative trait '", names(data)[j],
             "' must be numeric.")
      }
      
      if (any(!is.finite(x[!is.na(x)]))) {
        stop("Quantitative traits must contain only finite values or NA.")
      }
      
      valid <- !is.na(x)
      x_valid <- x[valid]
      g <- droplevels(clusters[valid])
      
      if (length(x_valid) < 2 || nlevels(g) < 2 ||
          length(unique(x_valid)) < 2) {
        scores[j] <- 0
        next
      }
      
      overall_mean <- mean(x_valid)
      total_ss <- sum((x_valid - overall_mean)^2)
      
      between_ss <- 0
      for (lev in levels(g)) {
        group_values <- x_valid[g == lev]
        between_ss <- between_ss +
          length(group_values) * (mean(group_values) - overall_mean)^2
      }
      
      scores[j] <- if (total_ss > 0) {
        min(1, max(0, between_ss / total_ss))
      } else {
        0
      }
      
    } else if (trait_type == "ordinal") {
      
      if (!is.ordered(x)) {
        stop("Ordinal trait '", names(data)[j],
             "' must be an ordered factor.")
      }
      
      ranks <- as.numeric(x)
      valid <- !is.na(ranks)
      rank_values <- ranks[valid]
      g <- droplevels(clusters[valid])
      
      if (length(rank_values) < 2 || nlevels(g) < 2 ||
          length(unique(rank_values)) < 2) {
        scores[j] <- 0
        next
      }
      
      overall_mean <- mean(rank_values)
      total_ss <- sum((rank_values - overall_mean)^2)
      
      between_ss <- 0
      for (lev in levels(g)) {
        group_values <- rank_values[g == lev]
        between_ss <- between_ss +
          length(group_values) * (mean(group_values) - overall_mean)^2
      }
      
      scores[j] <- if (total_ss > 0) {
        min(1, max(0, between_ss / total_ss))
      } else {
        0
      }
      
    } else {
      
      # Qualitative and binary traits
      valid <- !is.na(x)
      x_valid <- as.character(x[valid])
      g <- droplevels(clusters[valid])
      
      if (length(x_valid) < 2 || nlevels(g) < 2 ||
          length(unique(x_valid)) < 2) {
        scores[j] <- 0
        next
      }
      
      tab <- table(g, x_valid)
      
      if (nrow(tab) < 2 || ncol(tab) < 2) {
        scores[j] <- 0
        next
      }
      
      expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
      
      # Cramer's V squared: descriptive association of trait and cluster
      chi_squared <- sum((tab - expected)^2 / expected)
      denominator <- sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)
      
      scores[j] <- if (denominator > 0) {
        min(1, max(0, chi_squared / denominator))
      } else {
        0
      }
    }
  }
  
  result <- data.frame(
    trait = names(data),
    type = types,
    separation = scores,
    row.names = NULL
  )
  
  result <- result[order(result$separation, decreasing = TRUE,
                         na.last = TRUE), ]
  
  class(result) <- c("mixed_trait", "data.frame")
  result
}
