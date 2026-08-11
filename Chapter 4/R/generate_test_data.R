# generate_test_data.R ----------------------------------------------------
# Build a realistic-shaped test fixture for buffer_capture_analysis().
#
# Returns a list:
#   $inner   sf, single-row state polygon
#   $outer   sf, single-row state polygon  (NOT assumed to nest inner --
#                                           caller chooses overlapping or
#                                           adjacent states as desired)
#   $tracts  sf county polygons with a synthetic Pop column
#
# County populations are random (area * uniform density). If you need real
# numbers, swap in the tidycensus snippet at the bottom of the file.
#
# Note on nesting: real adjacent US states do not nest. The pipeline asks
# for nesting batches; pick states that nest (e.g. one inside the other --
# rare for real states), or wrap a small state with a buffer of the larger
# one, or just use it as-is if you only want to exercise the geometry on
# overlapping shapes. Three options are sketched below.
# -------------------------------------------------------------------------

generate_test_data <- function(
  state_inner   = "LA",
  state_outer   = "TX",
  year          = 2022,
  density_range = c(0.001, 0.005),
  seed          = 20250507
) {
  if (!requireNamespace("tigris", quietly = TRUE))
    stop("Install the 'tigris' package to use this generator.")

  options(tigris_use_cache = TRUE, tigris_progress = FALSE)

  states_raw   <- tigris::states(cb = TRUE, year = year)
  counties_raw <- tigris::counties(state = c(state_inner, state_outer),
                                   cb = TRUE, year = year)

  inner <- states_raw[states_raw$STUSPS == state_inner, ]
  outer <- states_raw[states_raw$STUSPS == state_outer, ]
  if (nrow(inner) == 0 || nrow(outer) == 0)
    stop("Could not find both states. Check state codes.")

  set.seed(seed)
  counties_proj <- sf::st_transform(counties_raw, 5070)   # CONUS Albers, m
  area_m2       <- as.numeric(sf::st_area(counties_proj))
  density       <- runif(nrow(counties_raw),
                         density_range[1], density_range[2])

  tracts <- counties_raw
  tracts$Pop <- round(area_m2 * density)

  list(inner = inner, outer = outer, tracts = tracts,
       state_inner = state_inner, state_outer = state_outer)
}

# Example: synthesise a strictly-nesting pair by buffering the inner state
# outward to make a fake "outer". Useful when you need real nesting for
# the pipeline's assumption.
make_nesting_outer <- function(inner, buffer_km = 200, proj_crs = 5070) {
  inner_p <- sf::st_transform(inner, proj_crs)
  buf     <- sf::st_buffer(inner_p, dist = buffer_km * 1000)
  buf$NAME <- paste0(inner$NAME[1], " + ", buffer_km, "km buffer")
  buf
}

# (Optional) real ACS population. Requires `tidycensus` and a Census API
# key. Replace `tracts$Pop` after running generate_test_data():
#
#   library(tidycensus)
#   pop <- tidycensus::get_acs(geography = "county", variables = "B01003_001",
#                              state = c("LA", "TX"), year = 2022,
#                              geometry = FALSE) |>
#     dplyr::select(GEOID, Pop = estimate)
#   td$tracts <- dplyr::left_join(td$tracts, pop, by = "GEOID")
