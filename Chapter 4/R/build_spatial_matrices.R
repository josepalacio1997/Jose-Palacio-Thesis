# build_spatial_matrices.R --------------------------------------------------
#
# Implements the full spatial-matrix construction of §2.2 of the paper
# "Inferring Spatial Transmission Dynamics of Respiratory Syncytial Virus
#  Across Houston Wastewater Treatment Plants" (Palacio et al., 2026).
#
# Given the WWTP polygons and a census-tract layer with population, returns
# the four matrices that appear in the paper:
#
#   D_tilde  — asymmetric Population Extended Hausdorff Distance
#              (paper eq. 3)
#   D        — symmetric PEHD, d_ij = max(d_tilde_ij, d_tilde_ji)
#              (paper eq. 4)
#   K        — Gaussian RBF kernel with zero diagonal
#              (paper eq. 5)
#   W        — row-normalized spatial weight matrix
#              (paper eq. 6)
#
# Everything is computed in one function call to keep the QMDs uncluttered
# and to make the pipeline reproducible with a single entry point.
#
# The PEHD algorithm (dot sampling, quantile-based effective distance, and
# symmetrization) is a direct implementation of the prototype code provided
# by J. Schedler; see also Schedler (2020) for the underlying Hausdorff-based
# framework. This file re-packages that prototype into a single function
# returning all four matrices used in the paper.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
})

#' Build the four spatial matrices (D_tilde, D, K, W) from WWTP polygons and
#' a census-tract layer with population.
#'
#' @param regions          sf polygon layer of WWTP service areas
#'                         (projected CRS in meters; e.g., UTM zone 15N).
#' @param tracts           sf polygon layer of census tracts (same CRS as
#'                         regions) with a numeric population column.
#' @param pop_col          Name of the population column in `tracts`.
#'                         Default: "POP24_5YR" (ACS 2020-2024 5-year).
#' @param people_per_dot   Monte Carlo resolution: each dot represents
#'                         ~ this many people (paper: nu = 200).
#' @param target_fraction  Population quantile that defines the PEHD.
#'                         Default 0.5 (population median).
#' @param sigma_K          Optional user-supplied bandwidth for the RBF
#'                         kernel (in km). If NULL (default), uses the
#'                         median-heuristic sigma^2 = median(d_ij^2) / 2.
#' @param labels           Character vector of length nrow(regions) used
#'                         as row/col names of the returned matrices.
#'                         Default: paste0("R", 1:m).
#'
#' @return A named list with:
#'   D_tilde         : m x m asymmetric PEHD matrix (km)
#'   D               : m x m symmetric PEHD matrix (km)
#'   K               : m x m Gaussian RBF kernel matrix (zero diagonal)
#'   W               : m x m row-normalized spatial weight matrix
#'   sigma_K         : bandwidth used (km)
#'   people_per_dot  : nu used
#'   target_fraction : quantile used (0.5 = median = PEHD)
#'   dots_per_region : list of length m with the sfc dot geometries
#'                     (useful for visualization; can be discarded).
build_spatial_matrices <- function(regions,
                                   tracts,
                                   pop_col         = "POP24_5YR",
                                   people_per_dot  = 200,
                                   target_fraction = 0.5,
                                   sigma_K         = NULL,
                                   labels          = NULL) {
  # -------- Input validation --------------------------------------------
  stopifnot(inherits(regions, "sf"), inherits(tracts, "sf"))
  if (isTRUE(sf::st_is_longlat(regions)))
    stop("`regions` is in a geographic CRS; project to a metric CRS first.")
  if (isTRUE(sf::st_is_longlat(tracts)))
    stop("`tracts` is in a geographic CRS; project to a metric CRS first.")
  if (sf::st_crs(tracts) != sf::st_crs(regions))
    stop("`tracts` and `regions` must share a CRS.")
  if (!pop_col %in% names(tracts))
    stop("Column '", pop_col, "' not found in tracts.")
  if (target_fraction <= 0 || target_fraction >= 1)
    stop("target_fraction must be in (0, 1).")
  if (people_per_dot <= 0)
    stop("people_per_dot must be positive.")

  m <- nrow(regions)
  if (is.null(labels)) labels <- paste0("R", seq_len(m))
  stopifnot(length(labels) == m)

  # =====================================================================
  # STEP 1 — Sample population dots per region (paper §2.2.1)
  # =====================================================================
  # For each region j, filter census tracts by centroid rule and sample
  # a regular grid of dots inside each tract with count round(P_c / nu).
  # Result: dots_per_region[[j]] is an sfc of ~ N_j / nu points, with
  # local density proportional to tract population.
  tract_centroids <- suppressWarnings(sf::st_centroid(tracts))

  dots_per_region <- vector("list", m)
  for (j in seq_len(m)) {
    keep <- lengths(sf::st_intersects(tract_centroids,
                                      regions[j, ])) > 0
    tj   <- tracts[keep, , drop = FALSE]

    pts_list <- lapply(seq_len(nrow(tj)), function(k) {
      pop <- tj[[pop_col]][k]
      if (is.na(pop) || pop <= 0) return(NULL)
      n_pts <- max(1L, as.integer(round(pop / people_per_dot)))
      tryCatch(sf::st_sample(tj[k, , drop = FALSE],
                             size = n_pts, type = "regular"),
               error = function(e) NULL)
    })
    parts <- Filter(function(x) !is.null(x) && length(x) > 0, pts_list)
    if (length(parts) == 0)
      stop("Region ", labels[j], " has no valid population dots.")
    dots_per_region[[j]] <- do.call(c, parts)
  }

  # =====================================================================
  # STEP 2 — Asymmetric PEHD matrix D_tilde (paper eq. 3)
  # =====================================================================
  # For each ordered pair (i, j), the distance from region i's polygon
  # to region j's dots. The population-median effective distance is the
  # k-th smallest, k = ceiling(target_fraction * n_dots). This is the
  # smallest r such that at least half of Omega_j's population lies
  # within distance r of Omega_i.
  D_tilde <- matrix(NA_real_, nrow = m, ncol = m,
                    dimnames = list(labels, labels))
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      dj    <- dots_per_region[[j]]
      dists <- as.numeric(sf::st_distance(regions[i, ], dj))
      n     <- length(dists)
      k     <- max(1L, as.integer(ceiling(target_fraction * n)))
      D_tilde[i, j] <- sort.int(dists, partial = k)[k]
    }
  }
  # Convert meters (sf default with UTM projection) to km.
  D_tilde <- D_tilde / 1000

  # =====================================================================
  # STEP 3 — Symmetric PEHD matrix D (paper eq. 4)
  # =====================================================================
  # d_ij = max(d_tilde_ij, d_tilde_ji). Guarantees symmetry, nonnegativity,
  # and zero diagonal — the three properties the downstream RBF + row
  # normalization require.
  D <- pmax(D_tilde, t(D_tilde))

  # =====================================================================
  # STEP 4 — Gaussian RBF kernel K (paper eq. 5)
  # =====================================================================
  # sigma_K^2 = 1/2 * median(d_ij^2 : i != j)  (Gretton et al. 2012)
  # K_ij = exp(-d_ij^2 / (2 sigma_K^2))
  if (is.null(sigma_K)) {
    off_sq  <- D[row(D) != col(D)]^2
    sigma_K <- sqrt(median(off_sq) / 2)
  }
  K <- exp(-D^2 / (2 * sigma_K^2))

  # =====================================================================
  # STEP 5 — Zero diagonal + row-normalize -> W (paper eq. 6)
  # =====================================================================
  # W_ii = 0 (self-contribution enters through (1 - rho) in the growth-rate
  # equation, so W must carry only neighbor contributions).
  # W_ij = K_ij / sum_k K_ik  — asymmetric because row sums differ across
  # plants (dense vs sparse neighborhoods).
  diag(K) <- 0
  rs      <- rowSums(K)
  if (any(rs == 0))
    stop("Row sum of K is zero for at least one region — pick a larger sigma_K.")
  W <- K / rs

  # -------- Sanity checks -----------------------------------------------
  stopifnot(all(abs(rowSums(W) - 1) < 1e-10),
            all(diag(W) == 0),
            isSymmetric(unname(D)),
            all(diag(D) == 0))

  # -------- Return -------------------------------------------------------
  list(
    D_tilde         = D_tilde,           # asymmetric PEHD (paper eq. 3)
    D               = D,                 # symmetric   PEHD (paper eq. 4)
    K               = K,                 # RBF kernel + zero diag (paper eq. 5)
    W               = W,                 # row-normalized W (paper eq. 6)
    sigma_K         = sigma_K,           # bandwidth used (km)
    people_per_dot  = people_per_dot,
    target_fraction = target_fraction,
    dots_per_region = dots_per_region    # for visualization / audits
  )
}
