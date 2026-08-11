# analysis.R --------------------------------------------------------------
# Top-level orchestrator. One call, one result list. Save helpers separate
# so callers control I/O.
# -------------------------------------------------------------------------

# Run the full buffer-capture analysis on a pair of nesting batches.
#
# Arguments
#   inner          sf, single-row polygon. The smaller batch; the buffer is
#                  centered on its centroid.
#   outer          sf, single-row polygon. The larger batch (assumed to
#                  contain `inner`). Population is sampled here.
#   tracts         sf polygon layer with a population column (defaults to
#                  "Pop"). Need only cover `outer` -- anything else is
#                  filtered out.
#   pop_col        Name of the population column in `tracts`.
#   people_per_dot Approx. people represented by each grid point. Smaller
#                  -> denser mesh, slower but more accurate.
#   target_fraction Fraction of dots the optimised buffer must capture.
#   tolerance      list(frac = ..., radius_m = ...). Bisection stopping
#                  criteria. Pass partial -- defaults fill in.
#   proj_crs       If non-NULL, all inputs are reprojected to this CRS
#                  before any geometry math. Use a projected CRS in
#                  meters (e.g. 5070 for CONUS, 32615 for Houston).
#   rule           Tract membership rule: "centroid_in" or "intersects".
#   colors, labels Plot styling -- partial overrides allowed (see plot.R
#                  for the full key set).
#   seed           For reproducibility of any tie-breaking inside spsample.
#
# Returns a list with: dots, centroid_guess, target, mesh, plots, inputs.
buffer_capture_analysis <- function(
  inner, outer, tracts,
  pop_col         = "Pop",
  people_per_dot  = 10,
  target_fraction = 0.95,
  tolerance       = list(frac = 0.005, radius_m = 25),
  proj_crs        = NULL,
  rule            = c("centroid_in", "intersects"),
  colors          = NULL,
  labels          = NULL,
  seed            = 20250507
) {
  rule      <- match.arg(rule)
  tolerance <- modifyList(list(frac = 0.005, radius_m = 25), tolerance %||% list())
  set.seed(seed)

  # 1. Reproject if requested ---------------------------------------------
  if (!is.null(proj_crs)) {
    inner  <- st_transform(inner, proj_crs)
    outer  <- st_transform(outer, proj_crs)
    tracts <- st_transform(tracts, proj_crs)
  }
  if (isTRUE(st_is_longlat(outer)))
    warning("Geometry is in a geographic CRS; buffer distances will be ",
            "interpreted as degrees, not meters. Set proj_crs.")

  # 2. Filter tracts inside outer -----------------------------------------
  tracts_in <- filter_tracts_in_batch(tracts, outer, rule = rule)
  if (nrow(tracts_in) == 0)
    stop("No tracts inside `outer` under rule '", rule, "'.")

  # 3. Population-proportional dots ---------------------------------------
  dots <- sample_population_dots(tracts_in,
                                 pop_col = pop_col,
                                 people_per_dot = people_per_dot)

  # 4. Centroid-radius guess ----------------------------------------------
  cguess <- centroid_radius_analysis(inner, outer, dots)

  # 5. Search for target-capture radius -----------------------------------
  ctarget <- find_capture_radius(inner, outer, dots,
                                 target_fraction = target_fraction,
                                 frac_tol        = tolerance$frac,
                                 radius_tol_m    = tolerance$radius_m)

  # 6. Mesh-density check -------------------------------------------------
  mesh <- mesh_density_report(outer, tracts_in, pop_col, people_per_dot,
                              target_radius = ctarget$radius)

  # 7. Plots --------------------------------------------------------------
  plots <- list(
    overview       = plot_overview(inner, outer, tracts_in,
                                   colors = colors, labels = labels),
    dots           = plot_dots(outer, inner, tracts_in, dots,
                               colors = colors, labels = labels),
    centroid_guess = plot_capture(
      inner, outer, dots, cguess$captured, cguess$buffer, cguess$intersection,
      cguess$center, cguess$cent_outer,
      title    = "Centroid-radius guess",
      subtitle = sprintf("r = %.2f km   |   captured = %.3f",
                         cguess$radius / 1000, cguess$proportion),
      colors = colors, labels = labels
    ),
    target_radius  = plot_capture(
      inner, outer, dots, ctarget$captured, ctarget$buffer, ctarget$intersection,
      ctarget$center, NULL,
      title    = sprintf("Radius for %.0f%% capture", 100 * target_fraction),
      subtitle = sprintf("r = %.2f km   |   achieved = %.3f",
                         ctarget$radius / 1000, ctarget$proportion),
      colors = colors, labels = labels
    )
  )

  list(
    inputs = list(
      inner = inner, outer = outer, tracts = tracts, tracts_in_outer = tracts_in,
      pop_col = pop_col, people_per_dot = people_per_dot,
      target_fraction = target_fraction, tolerance = tolerance,
      proj_crs = proj_crs, rule = rule
    ),
    dots           = dots,
    centroid_guess = cguess,
    target         = ctarget,
    mesh           = mesh,
    plots          = plots
  )
}

# Save figures and serialised result to disk (no-op if dir is NULL).
save_buffer_analysis <- function(result, fig_dir = NULL, out_dir = NULL,
                                 width = 9, height = 8, dpi = 200) {
  if (!is.null(fig_dir)) {
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    files <- c(overview       = "01_overview.png",
               dots           = "02_dots.png",
               centroid_guess = "03_centroid_guess.png",
               target_radius  = "04_target_radius.png")
    for (k in names(files)) {
      ggsave(file.path(fig_dir, files[[k]]), result$plots[[k]],
             width = width, height = height, dpi = dpi, bg = "white")
    }
    message("Wrote figures to ", fig_dir)
  }
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(result, file.path(out_dir, "result.rds"))
    write.csv(result$target$iter_log, file.path(out_dir, "iter_log.csv"),
              row.names = FALSE)
    message("Wrote outputs to ", out_dir)
  }
  invisible(result)
}

# Null-coalescing helper, in case it's not loaded elsewhere.
`%||%` <- function(a, b) if (!is.null(a)) a else b
