
#' Perform Principal Coordinates Analysis
#'
#' Performs Principal Coordinates Analysis (PCoA) using a
#' mixed-trait dissimilarity object.
#'
#' @param distance_object An object created by mixed_distance().
#' @param dimensions Number of coordinate dimensions to return.
#'
#' @return An object of class "mixed_pcoa" containing the
#'   coordinates, eigenvalues, and percentage of variation explained.
#'
#' @importFrom stats cmdscale
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
#' # Assign genotype names
#' rownames(rice_data) <- paste0("Rice_", sprintf("%02d", 1:20))
#'
#' # Specify trait types in column order
#' trait_types <- c(
#'   "quantitative", "quantitative", "qualitative",
#'   "binary", "ordinal"
#' )
#'
#' # Step 1: Calculate mixed-trait dissimilarities
#' dist_result <- mixed_distance(rice_data, trait_types)
#'
#' # Step 2: Perform Principal Coordinates Analysis
#' pcoa_result <- mixed_pcoa(
#'   dist_result,
#'   dimensions = 2
#' )
#'
#' # Step 3: Inspect the PCoA results
#' print(pcoa_result)
#'
#' # View the genotype coordinates
#' head(pcoa_result$coordinates)
#'
#' # View the percentage of variation explained
#' pcoa_result$percent_explained

#' @export
mixed_pcoa <- function(distance_object, dimensions = 2) {
  
  if (!inherits(distance_object, "mixed_distance")) {
    stop("'distance_object' must be created by mixed_distance().")
  }
  
  dmat <- distance_object$dissimilarity
  
  if (!is.matrix(dmat) || nrow(dmat) != ncol(dmat)) {
    stop("The dissimilarity matrix must be square.")
  }
  
  n <- nrow(dmat)
  
  if (n < 3) {
    stop("At least three observations are required for PCoA.")
  }
  
  if (anyNA(dmat) || any(!is.finite(dmat))) {
    stop("The dissimilarity matrix contains missing or non-finite values.")
  }
  
  if (any(dmat < 0)) {
    stop("Dissimilarities must be non-negative.")
  }
  
  if (!isTRUE(all.equal(unname(dmat), unname(t(dmat))))) {
    stop("The dissimilarity matrix must be symmetric.")
  }
  
  if (any(abs(diag(dmat)) > 1e-10)) {
    stop("The diagonal of the dissimilarity matrix must be zero.")
  }
  
  if (length(dimensions) != 1 ||
      !is.numeric(dimensions) ||
      is.na(dimensions) ||
      !is.finite(dimensions) ||
      dimensions != floor(dimensions) ||
      dimensions < 1 ||
      dimensions >= n) {
    stop("'dimensions' must be an integer from 1 to n - 1.")
  }
  
  pcoa_result <- stats::cmdscale(
    stats::as.dist(dmat),
    k = as.integer(dimensions),
    eig = TRUE,
    add = TRUE
  )
  
  coordinates <- pcoa_result$points
  eigenvalues <- pcoa_result$eig
  
  positive_eigenvalues <- pmax(eigenvalues, 0)
  total_positive <- sum(positive_eigenvalues)
  
  if (total_positive > 0) {
    percent_explained <- 100 * positive_eigenvalues / total_positive
  } else {
    percent_explained <- rep(NA_real_, length(eigenvalues))
  }
  
  result <- list(
    coordinates = coordinates,
    eigenvalues = eigenvalues,
    percent_explained = percent_explained,
    dimensions = as.integer(dimensions),
    correction = "Cailliez"
  )
  
  class(result) <- "mixed_pcoa"
  
  result
}
