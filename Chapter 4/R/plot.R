# plot.R ------------------------------------------------------------------
# All plotting helpers. They take sf objects + a colors list + a labels
# list and return ggplot objects. Saving is left to the caller.
#
# Color and label lists support partial overrides via modifyList -- pass
# only the keys you want to change.
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
})

# Default palette. Every key is documented so callers know what to override.
default_colors <- function() {
  list(
    inner_fill     = "#f5d4cf",   # inner batch fill
    inner_stroke   = "#a8331f",   # inner batch outline
    outer_fill     = "#eef5fb",   # outer batch fill
    outer_stroke   = "#1f6fa8",   # outer batch outline
    tract_stroke   = "grey75",    # tract outlines
    dot_outside    = "#1f6fa8",   # dot color when outside captured region
    dot_inside     = "#7a4a00",   # dot color when inside captured region
    buffer_stroke  = "#444444",   # buffer circle outline
    intersection   = "#f4cf6c",   # buffer ∩ outer fill
    centroid_inner = "#a8331f",   # inner centroid marker
    centroid_outer = "#1f6fa8",   # outer centroid marker
    centroid_line  = "grey30",    # line between the two centroids
    heatmap_low    = "#fff7ed",   # heatmap min (zero distance)
    heatmap_high   = "#9a3412",   # heatmap max
    heatmap_text   = "grey15"     # cell text
  )
}

default_labels <- function() {
  list(inner = "inner batch", outer = "outer batch")
}

# Merge user-supplied list with the defaults (user keys win).
.merge_with_defaults <- function(user, defaults) {
  if (is.null(user)) return(defaults)
  modifyList(defaults, user)
}

# ---- Plot 1: overview of the two batches with tracts inside outer ------
plot_overview <- function(inner, outer, tracts_in_outer,
                          colors = NULL, labels = NULL) {
  colors <- .merge_with_defaults(colors, default_colors())
  labels <- .merge_with_defaults(labels, default_labels())

  ggplot() +
    geom_sf(data = outer, fill = colors$outer_fill,
            color = colors$outer_stroke, linewidth = 0.6, alpha = 0.6) +
    geom_sf(data = inner, fill = colors$inner_fill,
            color = colors$inner_stroke, linewidth = 0.6, alpha = 0.6) +
    geom_sf(data = tracts_in_outer, fill = NA,
            color = colors$outer_stroke, linewidth = 0.25) +
    labs(title    = sprintf("%s and %s", labels$inner, labels$outer),
         subtitle = sprintf("Tracts shown for %s only", labels$outer),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11)
}

# ---- Plot 2: dots inside outer ----------------------------------------
plot_dots <- function(outer, inner, tracts, dots,
                      colors = NULL, labels = NULL) {
  colors <- .merge_with_defaults(colors, default_colors())
  labels <- .merge_with_defaults(labels, default_labels())

  ppd <- if ("tract_pop" %in% names(dots)) {
    # rough avg, in case people_per_dot wasn't stored
    round(sum(dots$tract_pop) / nrow(dots) / nrow(tracts))
  } else NA

  ggplot() +
    geom_sf(data = inner,  fill = colors$inner_fill,
            color = colors$outer_stroke, linewidth = 0.5) +
    geom_sf(data = outer,  fill = colors$outer_fill,
            color = colors$outer_stroke, linewidth = 0.5) +
    geom_sf(data = tracts, fill = NA, color = colors$tract_stroke,
            linewidth = 0.2) +
    geom_sf(data = dots, color = colors$outer_stroke,
            size = 0.18, alpha = 0.6) +
    labs(title    = sprintf("Population-proportional dots in %s", labels$outer),
         subtitle = sprintf("%d dots", nrow(dots)),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11)
}

# ---- Plot 3 / 4: buffer + intersection + captured/uncaptured dots ------
plot_capture <- function(inner, outer, dots, captured, buffer, intersection,
                         center_inner, center_outer = NULL,
                         title = "", subtitle = "",
                         colors = NULL, labels = NULL) {
  colors <- .merge_with_defaults(colors, default_colors())
  labels <- .merge_with_defaults(labels, default_labels())

  dots$captured <- captured

  p <- ggplot() +
    geom_sf(data = outer, fill = colors$outer_fill,
            color = colors$outer_stroke, linewidth = 0.5) +
    geom_sf(data = inner, fill = colors$inner_fill,
            color = colors$inner_stroke, linewidth = 0.5) +
    geom_sf(data = intersection, fill = colors$intersection,
            color = NA, alpha = 0.45) +
    geom_sf(data = buffer, fill = NA, color = colors$buffer_stroke,
            linewidth = 0.4, linetype = "22") +
    geom_sf(data = dots, aes(color = captured), size = 0.22, alpha = 0.75) +
    scale_color_manual(
      values = c(`FALSE` = colors$dot_outside, `TRUE` = colors$dot_inside),
      labels = c(`FALSE` = "outside", `TRUE` = "inside"),
      name   = "Dot status"
    ) +
    geom_sf(data = center_inner, color = colors$centroid_inner, size = 2.5)

  if (!is.null(center_outer)) {
    line <- st_sfc(st_linestring(rbind(
      st_coordinates(center_inner),
      st_coordinates(center_outer)
    )), crs = st_crs(inner))
    p <- p +
      geom_sf(data = line, color = colors$centroid_line, linewidth = 0.4) +
      geom_sf(data = center_outer, color = colors$centroid_outer, size = 2.5)
  }

  p +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# Heatmap of a pairwise distance matrix as returned by
# region_distance_matrix().
#
# Arguments
#   d            Numeric matrix with rownames/colnames (and ideally
#                `units` and `method` attributes).
#   show_values  Logical or NULL. If NULL (default), show numbers only
#                when the matrix has <= 20 rows (legible cell text).
#   digits       Decimal places for cell text.
#   cluster      If TRUE, reorder rows/cols by hierarchical clustering on
#                the distance matrix (helps see structure when N is big).
#   colors       Partial color override (see default_colors()).
#   title        Optional plot title; defaults to method + units.
plot_distance_heatmap <- function(d,
                                  show_values = NULL,
                                  digits      = 0,
                                  cluster     = FALSE,
                                  colors      = NULL,
                                  title       = NULL) {
  colors <- .merge_with_defaults(colors, default_colors())
  units  <- attr(d, "units")  %||% "km"
  method <- attr(d, "method") %||% "centroid"
  target_fraction <- attr(d, "target_fraction") 
  
  if (is.null(show_values)) show_values <- nrow(d) <= 20

  if (cluster && nrow(d) > 1) {
    ord <- stats::hclust(stats::as.dist(d))$order
    d   <- d[ord, ord, drop = FALSE]
  }

  df <- expand.grid(row = rownames(d), col = colnames(d),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$dist <- as.vector(d)
  df$row  <- factor(df$row, levels = rownames(d))
  df$col  <- factor(df$col, levels = colnames(d))

  p <- ggplot(df, aes(x = col, y = row, fill = dist)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low  = colors$heatmap_low,
                        high = colors$heatmap_high,
                        name = sprintf("Distance (%s)", units)) +
    scale_y_discrete(limits = rev) +
    coord_equal() +
    labs(title = title %||% .heatmap_default_title(method, units, target_fraction),
         x = if (method == "population") "target region (population)" else NULL,
         y = if (method == "population") "reference region (centroid)" else NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          panel.grid      = element_blank())

  if (show_values) {
    p <- p + geom_text(aes(label = formatC(dist, format = "f", digits = digits)),
                       size = 3, color = colors$heatmap_text)
  }
  p
}

# Build a default title that reflects the method.
.heatmap_default_title <- function(method, units, target_fraction = NA) {
  if (method == "population") {
    sprintf("Radius to capture %.0f%% of population (%s)",
            100 * target_fraction, units)
  } else {
    sprintf("Pairwise %s distance (%s)", method, units)
  }
}
