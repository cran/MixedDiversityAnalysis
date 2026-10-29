## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----load-package-------------------------------------------------------------
library(MixedDiversityAnalysis)

## ----example-data-------------------------------------------------------------
trait_data <- data.frame(
  Accession = paste0("G", 1:8),
  PlantHeight = c(85, 92, 78, 101, 88, 95, 82, 99),
  SeedWeight = c(24.5, 27.1, 22.8, 29.3, 25.6, 28.0, 23.9, 30.1),
  SeedColour = c(
    "Red", "Brown", "Red", "Black",
    "Brown", "Black", "Red", "Brown"
  ),
  DiseaseResistance = ordered(
    c(
      "High", "Medium", "Low", "High",
      "Medium", "High", "Low", "Medium"
    ),
    levels = c("Low", "Medium", "High")
  ),
  Lodging = c(
    "Absent", "Present", "Absent", "Present",
    "Absent", "Present", "Absent", "Present"
  ),
  stringsAsFactors = FALSE
)

trait_data

## ----prepare-data-------------------------------------------------------------
accession_ids <- trait_data$Accession

traits <- trait_data[, -1, drop = FALSE]
rownames(traits) <- accession_ids

trait_types <- c(
  "quantitative",
  "quantitative",
  "qualitative",
  "ordinal",
  "binary"
)

str(traits)
trait_types

## ----distance-analysis--------------------------------------------------------
dissimilarity <- mixed_distance(
  data = traits,
  types = trait_types
)

dissimilarity$dissimilarity

## ----hierarchical-clustering--------------------------------------------------
clustering <- mixed_cluster(
  distance_object = dissimilarity,
  k = 3,
  method = "average"
)

clustering$clusters
clustering$sizes
clustering$silhouette_mean

## ----pcoa-analysis------------------------------------------------------------
ordination <- mixed_pcoa(
  distance_object = dissimilarity,
  dimensions = 2
)

ordination$coordinates
ordination$percent_explained[1:2]

## ----evaluation-analysis------------------------------------------------------
cluster_results <- clustering$silhouette_mean

trait_results <- mixed_diversity(
  data = traits,
  types = trait_types,
  clusters = clustering$clusters
)

cluster_results
trait_results

## ----dendrogram-plot, fig.width=8, fig.height=6-------------------------------
mixed_visualize(
  object = clustering,
  type = "dendrogram",
  main = "Hierarchical Clustering of Crop Accessions"
)

## ----silhouette-plot, fig.width=8, fig.height=6-------------------------------
mixed_visualize(
  object = clustering,
  type = "silhouette",
  main = "Cluster Silhouette Analysis"
)

## ----pcoa-plot, fig.width=8, fig.height=6-------------------------------------
mixed_visualize(
  object = ordination,
  type = "pcoa",
  clusters = clustering$clusters,
  labels = accession_ids,
  main = "PCoA of Mixed-Trait Dissimilarities"
)

## ----trait-plot, fig.width=8, fig.height=6------------------------------------
mixed_visualize(
  object = trait_results,
  type = "trait",
  main = "Trait Separation Across Clusters"
)

