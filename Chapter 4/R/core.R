# core.R ------------------------------------------------------------------
# Pure analysis functions. No global state, no I/O. Each function takes sf
# objects and parameters in, returns sf objects or lists out.
#
# Geometric assumption: `inner` nests inside `outer`. Non-nesting cases
# should be handled by the caller before reaching here.
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(sp)
})

# Null-coalescing helper used across files.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Filter a tract layer to the rows that lie inside `batch`.
#   rule = "centroid_in" -> tract's centroid is inside batch (clean partition)
#   rule = "intersects"  -> tract overlaps batch at all (inclusive, may
#                           double-count tracts straddling a boundary if you
#                           call this twice with different batches)
filter_tracts_in_batch <- function(tracts, batch,
                                   rule = c("centroid_in", "intersects")) {
  rule <- match.arg(rule)
  if (rule == "centroid_in") {
    cents <- suppressWarnings(st_centroid(tracts))
    keep  <- lengths(st_intersects(cents, batch)) > 0
  } else {
    keep  <- lengths(st_intersects(tracts, batch)) > 0
  }
  tracts[keep, , drop = FALSE]
}

# Generate a population-proportional regular grid of points inside each
# tract. Each tract i gets round(Pop_i / people_per_dot) points via
# sp::spsample(type = "regular").
sample_population_dots <- function(tracts, pop_col = "Pop",
                                   people_per_dot = 10) {
  if (!pop_col %in% names(tracts))
    stop("Column '", pop_col, "' not found in tracts.")
  if (people_per_dot <= 0)
    stop("people_per_dot must be positive.")

  sample_one <- function(tract_row, idx) {
    pop <- tract_row[[pop_col]]
    if (is.na(pop) || pop <= 0) return(NULL)
    n_pts <- max(1L, as.integer(round(pop / people_per_dot)))
    poly_sp <- as(tract_row, "Spatial")
    pts <- tryCatch(sp::spsample(poly_sp, n = n_pts, type = "regular"),
                    error = function(e) NULL)
    if (is.null(pts) || length(pts) == 0)
      pts <- tryCatch(sp::spsample(poly_sp, n = n_pts + 2L, type = "regular"),
                      error = function(e) NULL)
    if (is.null(pts) || length(pts) == 0) return(NULL)
    out <- st_as_sf(pts)
    st_crs(out) <- st_crs(tracts)
    out$tract_idx <- idx
    out$tract_pop <- pop
    out
  }

  parts <- lapply(seq_len(nrow(tracts)), function(i)
    sample_one(tracts[i, , drop = FALSE], idx = i))
  out <- do.call(rbind, Filter(Negate(is.null), parts))
  if (is.null(out) || nrow(out) == 0)
    stop("No dots generated. Check pop_col / people_per_dot / tract geometry.")
  out
}

# Fraction of dots that fall inside (buffer(center, radius) ∩ outer).
# `center` is an sfc point, `dots` are sf points.
captured_fraction <- function(center, outer, dots, radius) {
  if (radius <= 0) return(0)
  buf <- st_buffer(center, dist = radius)
  inter <- st_intersection(buf, st_geometry(outer))
  if (length(inter) == 0) return(0)
  mean(lengths(st_intersects(dots, inter)) > 0)
}

# Centroid-radius guess: buffer centered on inner's centroid with radius
# equal to the distance to outer's centroid.
centroid_radius_analysis <- function(inner, outer, dots) {
  cent_inner <- suppressWarnings(st_centroid(st_geometry(inner)))
  cent_outer <- suppressWarnings(st_centroid(st_geometry(outer)))
  radius     <- as.numeric(st_distance(cent_inner, cent_outer))
  buf        <- st_sf(geometry = st_buffer(cent_inner, dist = radius))
  inter      <- st_intersection(buf, st_geometry(outer))
  hits       <- lengths(st_intersects(dots, inter)) > 0
  list(
    center       = cent_inner,
    cent_outer   = cent_outer,
    radius       = radius,
    buffer       = buf,
    intersection = inter,
    captured     = hits,
    proportion   = mean(hits)
  )
}

# Bisection for the radius (centered at inner's centroid) that captures
# `target_fraction` of the dots inside (buffer ∩ outer).
find_capture_radius <- function(inner, outer, dots,
                                target_fraction = 0.5,
                                frac_tol        = 0.005,
                                radius_tol_m    = 25,
                                max_iters       = 60) {
  cent_inner <- suppressWarnings(st_centroid(st_geometry(inner)))

  # Upper bound: distance to farthest vertex of outer, with 5% slack.
  outer_xy <- st_coordinates(st_cast(st_geometry(outer), "POINT"))
  cxy      <- st_coordinates(cent_inner)
  r_hi     <- 1.05 * sqrt(max((outer_xy[, "X"] - cxy[1])^2 +
                              (outer_xy[, "Y"] - cxy[2])^2))
  r_lo     <- 0

  cap_at <- function(r) captured_fraction(cent_inner, outer, dots, r)

  if (cap_at(r_hi) < target_fraction)
    stop("Cannot reach target_fraction=", target_fraction,
         " even at maximum radius. Check inputs.")

  log_rows <- vector("list", max_iters)
  iter <- 0L
  r_mid <- NA_real_; frac_mid <- NA_real_

  repeat {
    iter <- iter + 1L
    if (iter > max_iters) break
    r_mid <- 0.5 * (r_lo + r_hi)
    frac_mid <- cap_at(r_mid)
    log_rows[[iter]] <- data.frame(iter = iter, r_lo = r_lo, r_hi = r_hi,
                                   r_mid = r_mid, frac_mid = frac_mid)
    if (abs(frac_mid - target_fraction) < frac_tol &&
        (r_hi - r_lo) < radius_tol_m) break
    if (frac_mid < target_fraction) r_lo <- r_mid else r_hi <- r_mid
  }

  buf   <- st_sf(geometry = st_buffer(cent_inner, dist = r_mid))
  inter <- st_intersection(buf, st_geometry(outer))
  hits  <- lengths(st_intersects(dots, inter)) > 0

  list(
    target_fraction = target_fraction,
    radius          = r_mid,
    proportion      = frac_mid,
    iterations      = iter,
    iter_log        = do.call(rbind, log_rows[seq_len(iter)]),
    buffer          = buf,
    intersection    = inter,
    captured        = hits,
    center          = cent_inner
  )
}

# Mesh-density check. See README / qmd for the derivation.
#   eps_stat: tolerable stat SE on captured fraction
#   eps_bnd : tolerable boundary-discretisation error
mesh_density_report <- function(outer, tracts, pop_col, people_per_dot,
                                target_radius,
                                eps_stat = 0.005, eps_bnd = 0.005) {
  A_outer  <- as.numeric(sum(st_area(outer)))
  A_tracts <- as.numeric(sum(st_area(tracts)))
  total_pop <- sum(tracts[[pop_col]], na.rm = TRUE)
  density   <- total_pop / A_tracts
  s_current <- sqrt(people_per_dot / density)
  N_current <- total_pop / people_per_dot
  L         <- 2 * pi * target_radius

  SE_stat_current  <- 1 / (2 * sqrt(N_current))
  eps_bnd_current  <- (L * s_current) / (2 * A_tracts)

  k_max_stat <- total_pop * 4 * eps_stat^2
  s_max_bnd  <- 2 * A_tracts * eps_bnd / L
  k_max_bnd  <- s_max_bnd^2 * density
  k_min      <- min(k_max_stat, k_max_bnd)
  k_safe     <- k_min / 2

  list(
    geometry = list(A_outer_km2         = A_outer / 1e6,
                    A_tracts_km2        = A_tracts / 1e6,
                    total_pop           = total_pop,
                    avg_density_per_km2 = density * 1e6),
    current = list(people_per_dot = people_per_dot,
                   spacing_m      = s_current,
                   N_dots         = N_current,
                   SE_stat        = SE_stat_current,
                   eps_boundary   = eps_bnd_current),
    tolerances = list(eps_stat = eps_stat, eps_bnd = eps_bnd),
    recommended = list(k_max_for_stat = k_max_stat,
                       k_max_for_bnd  = k_max_bnd,
                       k_min          = k_min,
                       k_safe         = k_safe,
                       spacing_safe_m = sqrt(k_safe / density),
                       N_safe         = total_pop / k_safe),
    meets_target = people_per_dot <= k_min
  )
}

# Pairwise distance matrix between a set of regions.
#
# Arguments
#   regions    sf polygons (or points). Any number of rows.
#   label_col  Column name in `regions` to use as row/column labels. If
#              NULL, falls back to NAME, BatchName, or row index.
#   method     "centroid"   -> distance between centroids (default; fast,
#                               symmetric, the usual choice)
#              "boundary"   -> minimum distance between polygon boundaries
#                               (zero for overlapping/touching pairs)
#              "hausdorff"  -> Hausdorff distance (max-of-min boundary
#                               distance; sensitive to shape, not just
#                               location)
#              "population" -> radius around region_i's centroid that
#                               captures `target_fraction` of region_j's
#                               population. ASYMMETRIC; diagonal is
#                               non-zero and measures each region's own
#                               spread. Requires `tracts` and a population
#                               column.
#   units      "km" (default), "m", or "mi". CRS must be projected.
#
#   --- only used when method = "population" ---
#   tracts          sf polygon layer with a population column.
#   pop_col         Name of the population column. Default "Pop".
#   target_fraction Capture fraction in (0, 1). Default 0.95.
#   people_per_dot  Mesh density for the dot sampling. Default 10.
#   rule            Tract-in-region rule, as in buffer_capture_analysis.
#
# Returns a numeric matrix with `units`, `method`, `label_col`, and (for
# population) `target_fraction` attributes. Symmetric for geometric
# methods, asymmetric for "population".
region_distance_matrix <- function(
  regions,
  label_col = NULL,
  method    = c("centroid", "boundary", "hausdorff", "population"),
  units     = c("km", "m", "mi"),
  tracts          = NULL,
  pop_col         = "Pop",
  target_fraction = 0.95,
  people_per_dot  = 10,
  rule            = c("centroid_in", "intersects"),
  verbose         = FALSE
) {
  method <- match.arg(method)
  units  <- match.arg(units)
  rule   <- match.arg(rule)
  if (isTRUE(st_is_longlat(regions)))
    stop("`regions` is in a geographic CRS; project to a metric CRS first.")

  # Resolve labels
  labels <- if (!is.null(label_col)) {
    if (!label_col %in% names(regions))
      stop("Column '", label_col, "' not found in regions.")
    as.character(regions[[label_col]])
  } else if ("NAME" %in% names(regions)) {
    as.character(regions$NAME)
  } else if ("BatchName" %in% names(regions)) {
    as.character(regions$BatchName)
  } else {
    as.character(seq_len(nrow(regions)))
  }
  if (anyDuplicated(labels))
    labels <- make.unique(labels)

  conv <- switch(units, m = 1, km = 1/1000, mi = 1/1609.344)
  N    <- nrow(regions)

  if (method == "population") {
    if (is.null(tracts))
      stop("method = 'population' requires `tracts`.")
    if (isTRUE(st_is_longlat(tracts)))
      stop("`tracts` is in a geographic CRS; project to a metric CRS first.")
    if (st_crs(tracts) != st_crs(regions))
      stop("`tracts` and `regions` must share a CRS.")
    if (!pop_col %in% names(tracts))
      stop("Column '", pop_col, "' not found in tracts.")
    if (target_fraction <= 0 || target_fraction >= 1)
      stop("target_fraction must be in (0, 1).")

    # Sample dots inside each region once
    dots_per_region <- vector("list", N)
    pop_per_region  <- rep(NA_real_, N)
    for (j in seq_len(N)) {
      tj <- filter_tracts_in_batch(tracts, regions[j, , drop = FALSE], rule = rule)
      if (nrow(tj) == 0) {
        warning("Region ", labels[j], " contains no tracts; row/col will be NA.")
        next
      }
      pop_per_region[j] <- sum(tj[[pop_col]], na.rm = TRUE)
      dots_per_region[[j]] <- tryCatch(
        sample_population_dots(tj, pop_col = pop_col,
                               people_per_dot = people_per_dot),
        error = function(e) {
          warning("Failed to sample dots in region ", labels[j], ": ", conditionMessage(e))
          NULL
        }
      )
    }


    d <- matrix(NA_real_, nrow = N, ncol = N)
    for (i in seq_len(N)) {
      for (jj in seq_len(N)) {
        dj <- dots_per_region[[jj]]
        if (is.null(dj) || nrow(dj) == 0) next
        dists <- as.numeric(st_distance(regions[i, ], dj))
        # Smallest radius that captures at least target_fraction of dots:
        # the k-th smallest distance, k = ceiling(f * n). Matches the
        # bisection's "first r at which captured >= f" behaviour.
        n <- length(dists)
        k <- max(1L, as.integer(ceiling(target_fraction * n)))
        d[i, jj] <- sort.int(dists, partial = k)[k]
        if (verbose) {
          eff_ppd <- pop_per_region[jj] / nrow(dj)
          cat(sprintf("  %s -> %s : d = %.1f m | k = %d / %d dots | eff people/dot = %.1f\n",
                      labels[i], labels[jj],
                      d[i, jj], k, n, eff_ppd))
        }
      }
    }
    if (anyNA(d)) {
      bad <- labels[apply(is.na(d), 1, any) | apply(is.na(d), 2, any)]
      stop("Population distance is undefined for region(s): ",
           paste(unique(bad), collapse = ", "),
           ". Check that each region contains tracts with non-NA, positive ",
           pop_col, ". Re-run after fixing inputs.")
    }
    d <- pmax(d, t(d)) 
    d <- d * conv

  } else {
    geom <- if (method == "centroid") {
      suppressWarnings(st_centroid(st_geometry(regions)))
    } else {
      st_geometry(regions)
    }
    d <- if (method == "hausdorff") {
      st_distance(geom, geom, which = "Hausdorff")
    } else {
      st_distance(geom, geom)
    }
    d <- as.matrix(d) * conv
  }

  dimnames(d) <- list(labels, labels)
  attr(d, "units")     <- units
  attr(d, "method")    <- method
  attr(d, "label_col") <- label_col %||% "auto"
  if (method == "population")
    attr(d, "target_fraction") <- target_fraction
  d
}
