
#' Calculate Mixed-Trait Dissimilarity
#'
#' Calculates pairwise dissimilarities for datasets containing
#' quantitative, qualitative, binary, and ordinal traits.
#'
#' @param data A data frame with observations in rows and traits
#'   in columns.
#' @param types A character vector specifying the type of each
#'   trait: "quantitative", "qualitative", "binary", or "ordinal".
#'
#' @return An object of class "mixed_distance", containing:
#'   \describe{
#'     \item{dissimilarity}{A pairwise dissimilarity matrix.}
#'     \item{data}{The input data used in the calculation.}
#'     \item{types}{The specified type of each trait.}
#'     \item{method}{The name of the dissimilarity method.}
#'   }
#'
#' @details
#' This function implements a custom mixed-trait dissimilarity
#' method and does not use Gower dissimilarity.
#'
#' Quantitative differences are divided by the observed range.
#' Qualitative and binary traits receive a dissimilarity of zero
#' for matching values and one for different values. Binary traits
#' are treated symmetrically.
#'
#' Ordinal traits must be ordered factors. Differences between
#' their integer ranks are divided by the number of rank intervals.
#'
#' For each pair, the overall dissimilarity is the mean of the
#' available trait-specific dissimilarities. Missing comparisons
#' are omitted. If no trait can be compared for a pair, its
#' dissimilarity is NA. The diagonal is always zero.
#'
#' @examples
#' # Create an expanded dataset with mixed trait types
#' rice_data <- data.frame(
#'   PlantHeight = c(
#'     90, 100, 110, 95, 85,
#'     105, 115, 92, 98, 108,
#'     88, 102, 112, 97, 83,
#'     107, 118, 93, 101, 99
#'   ),
#'   GrainYield = c(
#'     25, 30, 35, 28, 22,
#'     32, 38, 26, 29, 34,
#'     24, 31, 36, 27, 21,
#'     33, 40, 25, 30, 28
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
#' # Assign genotype names
#' rownames(rice_data) <- paste0("Rice_", sprintf("%02d", 1:20))
#'
#' # Specify the type of each trait in column order
#' trait_types <- c(
#'   "quantitative",
#'   "quantitative",
#'   "qualitative",
#'   "binary",
#'   "ordinal"
#' )
#'
#' # Calculate mixed-trait dissimilarities
#' result <- mixed_distance(rice_data, trait_types)
#'
#' # Display the dissimilarity matrix
#' round(result$dissimilarity, 3)

#' @export
mixed_distance <- function(data, types) {
  
  # 1. Validate the input data
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  
  if (ncol(data) == 0) {
    stop("'data' must contain at least one trait.")
  }
  
  if (nrow(data) < 2) {
    stop("'data' must contain at least two observations.")
  }
  
  if (!is.character(types)) {
    stop("'types' must be a character vector.")
  }
  
  if (length(types) != ncol(data)) {
    stop("'types' must have one entry for each column.")
  }
  
  if (anyNA(types)) {
    stop("'types' must not contain missing values.")
  }
  
  allowed_types <- c(
    "quantitative",
    "qualitative",
    "binary",
    "ordinal"
  )
  
  if (!all(types %in% allowed_types)) {
    stop(
      "Trait types must be quantitative, qualitative, binary, or ordinal."
    )
  }
  
  # Work on a copy so the original data is not modified
  prepared_data <- data
  
  # 2. Validate and prepare each trait
  for (k in seq_along(types)) {
    
    x <- prepared_data[[k]]
    observed <- x[!is.na(x)]
    
    if (types[k] == "quantitative") {
      
      if (!is.numeric(x) || is.factor(x)) {
        stop("Quantitative traits must be numeric.")
      }
      
      if (length(observed) == 0) {
        stop(
          "Quantitative traits must not contain only missing values."
        )
      }
      
      if (any(!is.finite(observed))) {
        stop(
          "Quantitative traits must not contain Inf or -Inf."
        )
      }
      
    } else if (types[k] == "qualitative") {
      
      if (!(is.factor(x) || is.character(x))) {
        stop(
          "Qualitative traits must be factors or character vectors."
        )
      }
      
      prepared_data[[k]] <- as.character(x)
      
    } else if (types[k] == "binary") {
      
      if (!(is.factor(x) || is.character(x) ||
            is.logical(x) || is.numeric(x))) {
        stop(
          "Binary traits must be factor, character, logical, or numeric."
        )
      }
      
      if (is.numeric(x) &&
          any(!is.finite(observed))) {
        stop(
          "Numeric binary traits must not contain Inf or -Inf."
        )
      }
      
      if (length(unique(observed)) != 2) {
        stop(
          "Each binary trait must have exactly two observed values."
        )
      }
      
    } else if (types[k] == "ordinal") {
      
      if (!is.ordered(x)) {
        stop(
          "Ordinal traits must be supplied as ordered factors."
        )
      }
      
      if (length(observed) == 0) {
        stop(
          "Ordinal traits must not contain only missing values."
        )
      }
      
      if (nlevels(x) < 2) {
        stop(
          "Ordinal traits must have at least two ordered levels."
        )
      }
    }
  }
  
  # 3. Calculate quantitative trait ranges
  trait_ranges <- rep(NA_real_, ncol(prepared_data))
  
  for (k in seq_along(types)) {
    if (types[k] == "quantitative") {
      trait_ranges[k] <- diff(
        range(prepared_data[[k]], na.rm = TRUE)
      )
    }
  }
  
  # 4. Initialize the dissimilarity matrix
  n <- nrow(prepared_data)
  p <- ncol(prepared_data)
  
  distance_matrix <- matrix(
    NA_real_,
    nrow = n,
    ncol = n
  )
  
  # 5. Calculate pairwise dissimilarities
  for (i in seq_len(n)) {
    for (j in i:n) {
      
      # An observation has zero distance from itself
      if (i == j) {
        distance_matrix[i, j] <- 0
        next
      }
      
      trait_distances <- rep(NA_real_, p)
      
      for (k in seq_len(p)) {
        
        x <- prepared_data[[k]][i]
        y <- prepared_data[[k]][j]
        
        # Skip comparisons involving missing values
        if (is.na(x) || is.na(y)) {
          next
        }
        
        if (types[k] == "quantitative") {
          
          r <- trait_ranges[k]
          
          if (r == 0) {
            trait_distances[k] <- 0
          } else {
            trait_distances[k] <- abs(x - y) / r
          }
          
        } else if (types[k] == "ordinal") {
          
          ranks <- as.integer(prepared_data[[k]])
          rank_range <- nlevels(prepared_data[[k]]) - 1
          
          trait_distances[k] <-
            abs(ranks[i] - ranks[j]) / rank_range
          
        } else {
          
          # Qualitative and symmetric binary traits
          trait_distances[k] <- as.numeric(x != y)
          
        }
      }
      
      valid <- !is.na(trait_distances)
      
      pair_distance <- if (any(valid)) {
        mean(trait_distances[valid])
      } else {
        NA_real_
      }
      
      # Store both directions
      distance_matrix[i, j] <- pair_distance
      distance_matrix[j, i] <- pair_distance
    }
  }
  
  # 6. Assign observation names
  observation_names <- rownames(prepared_data)
  
  if (is.null(observation_names)) {
    observation_names <- as.character(seq_len(n))
  }
  
  rownames(distance_matrix) <- observation_names
  colnames(distance_matrix) <- observation_names
  
  # 7. Return a structured object
  result <- list(
    dissimilarity = distance_matrix,
    data = prepared_data,
    types = types,
    method = "custom_mixed_dissimilarity"
  )
  
  class(result) <- "mixed_distance"
  
  return(result)
}




