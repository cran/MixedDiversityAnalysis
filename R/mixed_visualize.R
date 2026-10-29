
#' Visualize Mixed-Trait Clustering Results
#'
#' Creates graphical summaries of hierarchical clustering,
#' silhouette widths, principal coordinates analysis (PCoA),
#' and trait contributions to cluster separation.
#'
#' @param object An object returned by `mixed_cluster()`,
#'   `mixed_pcoa()`, or `mixed_diversity()`, depending on `type`.
#' @param type Character string specifying the plot type:
#'   `"dendrogram"`, `"silhouette"`, `"pcoa"`, or `"trait"`.
#' @param clusters Optional vector assigning observations to clusters.
#'   Used to colour PCoA points and optionally draw group ellipses.
#' @param labels Optional vector of observation labels for the PCoA plot.
#'   Must contain one non-missing value per observation.
#' @param main Optional character string specifying the plot title.
#' @param ellipse Logical; if `TRUE`, attempts to draw 95% group
#'   ellipses on the PCoA plot when cluster assignments are supplied.
#'   Ellipses are omitted for groups with insufficient or degenerate
#'   coordinate data.
#' @param ... Additional graphical parameters passed to the relevant
#'   base plotting function.
#'
#' @return Invisibly returns the supplied `object`. The function is
#'   called primarily for its graphical side effect.
#'
#' @details
#' The `"dendrogram"` option displays the hierarchical clustering tree
#' and highlights the resulting clusters. The `"silhouette"` option
#' displays observation-level silhouette widths and the mean silhouette
#' width. The `"pcoa"` option plots the first two principal coordinate
#' dimensions and can colour observations by cluster. Optional ellipses
#' summarize the spread of groups in the plotted coordinate space.
#' The `"trait"` option displays trait separation scores calculated by
#' `mixed_diversity()`.
#'
#' Plot appearance uses base R graphics. The function restores the
#' previous graphical parameters on exit.
#'
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
#' # Specify trait types
#' trait_types <- c(
#'   "quantitative", "quantitative", "qualitative",
#'   "binary", "ordinal"
#' )
#'
#' # Calculate mixed-trait dissimilarities
#' dist_result <- mixed_distance(rice_data, trait_types)
#'
#' # Perform hierarchical clustering
#' cluster_result <- mixed_cluster(
#'   dist_result,
#'   k = 3,
#'   method = "average"
#' )
#'
#' # Plot the dendrogram
#' mixed_visualize(
#'   object = cluster_result,
#'   type = "dendrogram",
#'   main = "Rice Genotype Clustering"
#' )
#'
#' # Plot the PCoA ordination
#' pcoa_result <- mixed_pcoa(dist_result, dimensions = 2)
#' mixed_visualize(
#'   object = pcoa_result,
#'   type = "pcoa",
#'   clusters = cluster_result$clusters,
#'   labels = rownames(rice_data),
#'   main = "PCoA of Rice Genotypes"
#' )
#'
#' # Plot trait separation scores
#' trait_result <- mixed_diversity(
#'   data = rice_data,
#'   types = trait_types,
#'   clusters = cluster_result$clusters
#' )
#' mixed_visualize(
#'   object = trait_result,
#'   type = "trait",
#'   main = "Trait Contributions to Cluster Separation"
#' )

#' @export
mixed_visualize <- function(
    object,
    type = c("dendrogram", "silhouette", "pcoa", "trait"),
    clusters = NULL,
    labels = NULL,
    main = NULL,
    ellipse = FALSE,
    ...
) {
  
  type <- match.arg(type)
  
  if (!is.logical(ellipse) || length(ellipse) != 1L ||
      is.na(ellipse)) {
    stop("'ellipse' must be TRUE or FALSE.")
  }
  
  palette <- c(
    "#2878B5", "#E07A35", "#26966F", "#8A63B8",
    "#C94C4C", "#7A8793", "#B59B25", "#D66BA0"
  )
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::par(
    bg = "white",
    fg = "#263238",
    col.axis = "#37474F",
    col.lab = "#263238",
    col.main = "#172B4D",
    font.main = 2,
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.85,
    las = 1
  )
  
  # ==================================================
  # 1. DENDROGRAM
  # ==================================================
  if (type == "dendrogram") {
    
    if (!inherits(object, "mixed_cluster")) {
      stop("'object' must be returned by mixed_cluster().")
    }
    
    graphics::par(mar = c(5, 5, 5, 2) + 0.1)
    
    graphics::plot(
      object$tree,
      hang = -1,
      main = if (is.null(main)) {
        "Hierarchical Clustering"
      } else {
        main
      },
      sub = "",
      xlab = "Observations",
      ylab = "Dissimilarity",
      cex = 0.9,
      col = "#455A64",
      ...
    )
    
    stats::rect.hclust(
      object$tree,
      k = object$k,
      border = palette[
        seq_len(min(object$k, length(palette)))
      ]
    )
    
    graphics::mtext(
      paste0(
        "Clusters: ", object$k,
        " | Linkage: ", object$method,
        " | Mean silhouette: ",
        sprintf("%.3f", object$silhouette_mean)
      ),
      side = 3,
      line = 0.35,
      cex = 0.7,
      col = "#546E7A"
    )
    
    # ==================================================
    # 2. SILHOUETTE PLOT
    # ==================================================
  } else if (type == "silhouette") {
    
    if (!inherits(object, "mixed_cluster")) {
      stop("'object' must be returned by mixed_cluster().")
    }
    
    sil <- object$silhouette
    widths <- sil[, "sil_width"]
    groups <- sil[, "cluster"]
    
    obs_labels <- rownames(sil)
    
    if (is.null(obs_labels)) {
      obs_labels <- as.character(seq_along(widths))
    }
    
    ord <- order(groups, -widths)
    
    widths <- widths[ord]
    groups <- groups[ord]
    obs_labels <- obs_labels[ord]
    
    group_levels <- sort(unique(groups))
    group_index <- match(groups, group_levels)
    
    bar_col <- palette[
      ((group_index - 1) %% length(palette)) + 1
    ]
    
    graphics::par(mar = c(5, 5, 4, 2) + 0.1)
    
    n_obs <- length(widths)
    ypos <- seq_len(n_obs)
    
    graphics::plot(
      NA,
      xlim = c(min(-0.08, min(widths) - 0.04), 1),
      ylim = c(0.5, n_obs + 0.5),
      xlab = "Silhouette width",
      ylab = "Observations",
      yaxt = "n",
      main = if (is.null(main)) {
        "Silhouette Analysis"
      } else {
        main
      },
      ...
    )
    
    graphics::abline(
      v = 0,
      col = "#9AA5B1",
      lty = 2
    )
    
    for (i in seq_along(widths)) {
      graphics::rect(
        xleft = 0,
        ybottom = ypos[i] - 0.38,
        xright = widths[i],
        ytop = ypos[i] + 0.38,
        col = bar_col[i],
        border = NA
      )
    }
    
    graphics::axis(
      side = 2,
      at = ypos,
      labels = obs_labels,
      tick = FALSE,
      cex.axis = if (n_obs > 25) 0.65 else 0.85
    )
    
    graphics::abline(
      v = object$silhouette_mean,
      col = "#B22222",
      lwd = 2,
      lty = 2
    )
    
    graphics::legend(
      "bottomright",
      legend = paste0(
        "Mean silhouette = ",
        sprintf("%.3f", object$silhouette_mean)
      ),
      col = "#B22222",
      lty = 2,
      lwd = 2,
      bty = "n",
      cex = 0.8
    )
    
    # ==================================================
    # 3. PCoA PLOT WITH OPTIONAL GROUP ELLIPSES
    # ==================================================
  } else if (type == "pcoa") {
    
    if (!inherits(object, "mixed_pcoa")) {
      stop("'object' must be returned by mixed_pcoa().")
    }
    
    coordinates <- as.matrix(object$coordinates)
    
    if (ncol(coordinates) < 2) {
      stop("At least two PCoA dimensions are required.")
    }
    
    n_obs <- nrow(coordinates)
    
    if (n_obs < 2 || any(!is.finite(coordinates[, 1:2]))) {
      stop("The first two PCoA coordinates must be finite.")
    }
    
    x <- coordinates[, 1]
    y <- coordinates[, 2]
    
    if (is.null(labels)) {
      labels <- rownames(coordinates)
      
      if (is.null(labels)) {
        labels <- as.character(seq_len(n_obs))
      }
    }
    
    if (length(labels) != n_obs || anyNA(labels)) {
      stop("'labels' must contain one non-missing value per observation.")
    }
    
    # Set up groups
    if (is.null(clusters)) {
      
      group <- factor(rep("Observations", n_obs))
      
    } else {
      
      if (length(clusters) != n_obs || anyNA(clusters)) {
        stop(
          "'clusters' must have one non-missing value per observation."
        )
      }
      
      group <- droplevels(as.factor(clusters))
      
      if (any(table(group) == 0)) {
        stop("Each cluster must contain at least one observation.")
      }
    }
    
    group_levels <- levels(group)
    
    group_index <- match(group, group_levels)
    
    point_col <- palette[
      ((group_index - 1) %% length(palette)) + 1
    ]
    
    # Axis labels
    pct <- object$percent_explained
    
    axis1 <- "PCoA 1"
    axis2 <- "PCoA 2"
    
    if (length(pct) >= 2 && all(is.finite(pct[1:2]))) {
      axis1 <- sprintf("PCoA 1 (%.2f%%)", pct[1])
      axis2 <- sprintf("PCoA 2 (%.2f%%)", pct[2])
    }
    
    # Helper to provide padding around plotted values
    expand_range <- function(z) {
      
      r <- range(z, finite = TRUE)
      span <- diff(r)
      
      pad <- if (span == 0) {
        max(abs(r), 1) * 0.05
      } else {
        span * 0.12
      }
      
      r + c(-pad, pad)
    }
    
    
    # Calculate each ellipse once and reuse it for plotting.
    # At least three distinct, non-collinear points are required.
    ellipse_list <- vector("list", length(group_levels))
    names(ellipse_list) <- group_levels
    
    if (ellipse && !is.null(clusters)) {
      
      for (g in group_levels) {
        
        idx <- which(group == g)
        
        xy <- cbind(x[idx], y[idx])
        xy <- unique(xy)
        
        # Skip groups with fewer than three unique points.
        if (nrow(xy) < 3) {
          next
        }
        
        # Calculate covariance safely.
        cov_matrix <- tryCatch(
          stats::cov(xy),
          error = function(e) NULL
        )
        
        if (is.null(cov_matrix) ||
            any(!is.finite(cov_matrix))) {
          next
        }
        
        # Calculate eigenvalues and eigenvectors safely.
        eig <- tryCatch(
          eigen(cov_matrix, symmetric = TRUE),
          error = function(e) NULL
        )
        
        if (is.null(eig) ||
            any(!is.finite(eig$values)) ||
            min(eig$values) <= 1e-12) {
          next
        }
        
        # Generate points for a 95% confidence ellipse.
        angles <- seq(0, 2 * pi, length.out = 200)
        
        radius <- sqrt(stats::qchisq(0.95, df = 2))
        
        ellipse_xy <- t(
          eig$vectors %*%
            diag(sqrt(eig$values) * radius) %*%
            rbind(cos(angles), sin(angles))
        )
        
        # Centre the ellipse using the same unique points
        # used to calculate the covariance matrix.
        ellipse_xy[, 1] <- ellipse_xy[, 1] + mean(xy[, 1])
        ellipse_xy[, 2] <- ellipse_xy[, 2] + mean(xy[, 2])
        
        # Store the ellipse for later plotting.
        ellipse_list[[g]] <- ellipse_xy
      }
    }
    
    
    # Include ellipse coordinates in plot limits
    xlim <- expand_range(x)
    ylim <- expand_range(y)
    
    valid_ellipses <- ellipse_list[
      !vapply(ellipse_list, is.null, logical(1))
    ]
    
    if (length(valid_ellipses) > 0) {
      
      ellipse_x <- unlist(lapply(
        valid_ellipses, function(z) z[, 1]
      ))
      
      ellipse_y <- unlist(lapply(
        valid_ellipses, function(z) z[, 2]
      ))
      
      xlim <- range(xlim, ellipse_x, finite = TRUE)
      ylim <- range(ylim, ellipse_y, finite = TRUE)
    }
    
    graphics::par(mar = c(5, 5, 5, 2) + 0.1)
    
    graphics::plot(
      x,
      y,
      xlim = xlim,
      ylim = ylim,
      type = "n",
      xlab = axis1,
      ylab = axis2,
      main = if (is.null(main)) {
        "Principal Coordinates Analysis"
      } else {
        main
      },
      ...
    )
    
    graphics::abline(
      h = 0,
      v = 0,
      col = "#E0E5EA",
      lty = 2
    )
    
    # Draw group ellipses before the points
    if (length(valid_ellipses) > 0) {
      
      for (g in names(valid_ellipses)) {
        
        ellipse_xy <- valid_ellipses[[g]]
        idx_g <- match(g, group_levels)
        
        graphics::lines(
          ellipse_xy[, 1],
          ellipse_xy[, 2],
          col = palette[
            ((idx_g - 1) %% length(palette)) + 1
          ],
          lwd = 1.5
        )
      }
    }
    
    # Draw observations
    graphics::points(
      x,
      y,
      pch = 21,
      bg = point_col,
      col = "white",
      cex = 1.5
    )
    
    graphics::text(
      x,
      y,
      labels = labels,
      pos = 3,
      offset = 0.5,
      cex = 0.85
    )
    
    # Cluster legend
    if (!is.null(clusters)) {
      
      graphics::legend(
        "topright",
        legend = group_levels,
        pt.bg = palette[
          ((seq_along(group_levels) - 1) %% length(palette)) + 1
        ],
        pch = 21,
        col = "white",
        bty = "n",
        title = "Cluster",
        cex = 0.8
      )
    }
    
    # Report ellipse availability
    if (ellipse && !is.null(clusters) &&
        length(valid_ellipses) < length(group_levels)) {
      
      graphics::mtext(
        "Some ellipses omitted: insufficient or degenerate group coordinates",
        side = 1,
        line = 3.5,
        cex = 0.65,
        col = "#7A5B00"
      )
    }
    
    # Axis percentages
    if (length(pct) >= 2 && all(is.finite(pct[1:2]))) {
      
      graphics::mtext(
        sprintf(
          "PCoA 1: %.2f%% | PCoA 2: %.2f%%",
          pct[1],
          pct[2]
        ),
        side = 3,
        line = 0.3,
        cex = 0.75,
        col = "#546E7A"
      )
    }
    
    # ==================================================
    # 4. TRAIT SEPARATION PLOT
    # ==================================================
  } else if (type == "trait") {
    
    if (!inherits(object, "mixed_trait")) {
      stop("'object' must be returned by mixed_diversity().")
    }
    
    required <- c("trait", "type", "separation")
    
    if (!all(required %in% names(object))) {
      stop(
        "'object' must contain trait, type, and separation columns."
      )
    }
    
    plot_data <- object[
      !is.na(object$separation),
      ,
      drop = FALSE
    ]
    
    if (nrow(plot_data) == 0) {
      stop("No non-missing trait separation scores are available.")
    }
    
    plot_data <- plot_data[
      order(plot_data$separation),
      ,
      drop = FALSE
    ]
    
    bar_colors <- ifelse(
      plot_data$type == "quantitative", "#2878B5",
      ifelse(
        plot_data$type == "qualitative", "#26966F",
        ifelse(
          plot_data$type == "binary", "#E07A35", "#8A63B8"
        )
      )
    )
    
    graphics::par(mar = c(7, 8, 5, 3) + 0.1)
    
    positions <- seq_len(nrow(plot_data))
    
    plot_title <- if (is.null(main)) {
      "Trait Contributions to Cluster Separation"
    } else {
      main
    }
    
    graphics::plot(
      NA,
      xlim = c(0, 1.12),
      ylim = c(0.5, nrow(plot_data) + 0.5),
      xlab = "Cluster separation score",
      ylab = "",
      yaxt = "n",
      xaxt = "n",
      main = plot_title,
      xaxs = "i",
      yaxs = "i"
    )
    
    graphics::axis(
      side = 1,
      at = seq(0, 1, by = 0.2),
      labels = sprintf("%.1f", seq(0, 1, by = 0.2))
    )
    
    graphics::axis(
      side = 2,
      at = positions,
      labels = plot_data$trait,
      las = 1,
      tick = FALSE
    )
    
    graphics::abline(
      v = seq(0.2, 1, by = 0.2),
      col = "#E6EAEE",
      lty = 3
    )
    
    graphics::rect(
      xleft = 0,
      ybottom = positions - 0.32,
      xright = plot_data$separation,
      ytop = positions + 0.32,
      col = bar_colors,
      border = NA
    )
    
    graphics::text(
      x = plot_data$separation + 0.025,
      y = positions,
      labels = sprintf("%.3f", plot_data$separation),
      adj = c(0, 0.5),
      cex = 0.9,
      col = "#263238"
    )
    
    graphics::legend(
      "bottom",
      inset = c(0, -0.30),
      xpd = NA,
      horiz = TRUE,
      legend = c(
        "Quantitative",
        "Qualitative",
        "Binary",
        "Ordinal"
      ),
      fill = c(
        "#2878B5",
        "#26966F",
        "#E07A35",
        "#8A63B8"
      ),
      border = NA,
      bty = "n",
      cex = 0.7,
      title = "Trait type"
    )
  }
  
  invisible(object)
}
