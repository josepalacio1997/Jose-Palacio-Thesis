# ============================================================================
# run_A_fsv_est_HD.R
#
# HAUSDORFF ABLATION variant of run_A_fsv_est.R
# -----------------------------------------------------------------------------
# Identical to run_A_fsv_est.R except the spatial weight matrix W is loaded
# from W_HD.rds (classical (sup-based) Hausdorff distance + RBF kernel +
# row-normalization), in place of the PEHD-derived W_pop50 used in the
# primary fit. All other model components, priors, HMC settings, and Stan
# files remain unchanged.
#
# Prerequisite: run paper_pipeline/01_W_matrix/compute_W_matrices.R first
# to generate ../results/W_HD.rds (and W_PEHD.rds for reference).
#
# In-sample fit del modelo FSV (freesigma_vec) con RHO ESTIMADO.
# Llamada: Rscript run_A_fsv_est_HD.R [ALPHA_RHO] [BETA_RHO] [WARMUP] [SAMPLING]
#   ej. Rscript run_A_fsv_est_HD.R                  # defaults: Beta(2,2)
#   ej. Rscript run_A_fsv_est_HD.R 2 2 4000 6000    # explicit
#
# Stan:   stan/renewal_A_fsv_est.stan       (unchanged)
# Output: outputs/stan_fits_fsv_rhoest_HD/rho_estimated.rds
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
ALPHA_RHO    <- if (length(args) >= 1) as.numeric(args[1]) else 2.0
BETA_RHO     <- if (length(args) >= 2) as.numeric(args[2]) else 2.0
arg_warmup   <- if (length(args) >= 3) suppressWarnings(as.integer(args[3])) else NA_integer_
arg_sampling <- if (length(args) >= 4) suppressWarnings(as.integer(args[4])) else NA_integer_
stopifnot(ALPHA_RHO > 0, BETA_RHO > 0)
JOB_NAME <- "rho_estimated"
cat(sprintf("\n=== IN-SAMPLE FSV (rho ESTIMATED, Beta(%g, %g)) ===\n",
            ALPHA_RHO, BETA_RHO))

suppressPackageStartupMessages({
  library(tidyverse); library(lubridate); library(sf)
  library(cmdstanr); library(posterior); library(purrr)
})

set.seed(42)
proj_root <- normalizePath(".", mustWork = FALSE)
cat("proj_root =", proj_root, "\n")

# ---- Config -----------------------------------------------------------------
n_chains      <- 3
iter_warmup   <- if (!is.na(arg_warmup))   arg_warmup   else 4000
iter_sampling <- if (!is.na(arg_sampling)) arg_sampling else 6000
adapt_delta   <- 0.95
max_treedepth <- 15
cat(sprintf("Config: %d chains x %d warmup x %d sampling (total=%d)\n",
            n_chains, iter_warmup, iter_sampling, iter_warmup + iter_sampling))

Batch_focus <- paste0("WWTP", 1:15)
raw_csv_path <- file.path(proj_root, "data",
                          "all_ts_observed_ILI_20250526_pool.csv")
shp_batch    <- file.path(proj_root, "data", "WWTP_Batch",
                          "WWTP_061621_mergedWWTPs to Batch_updatedon07192024.shp")
csv_ct_pop   <- file.path(proj_root, "data", "Census_Tracts",
                          "WWTP Batch Service Areas_ Census Tracts_Population.csv")
tracts_shp_path <- file.path(proj_root, "data", "Census_Tracts_Population",
                             "WWTP Batch Service Areas_Census Tracts_Population.shp")

casos_houston_semana <- 324
burnin               <- 15
cutoff_date_min      <- as.Date("2023-01-02")
prior_mean_gen_days  <- 7.5; prior_sd_gen_days <- 2.1
mean_shed_days       <- 4.6; sd_shed_days      <- 2.0

source(file.path(proj_root, "R", "core.R"))
source(file.path(proj_root, "R", "plot.R"))
source(file.path(proj_root, "R", "analysis.R"))

to_gamma <- function(mean_days, sd_days) {
  m <- mean_days / 7; s <- sd_days / 7
  list(shape = (m^2) / (s^2), scale = (s^2) / m)
}
make_shedding_kernel <- function(mean_shed_days, sd_shed_days, K = 2) {
  par_d   <- to_gamma(mean_shed_days, sd_shed_days)
  edges_d <- c(-0.5, seq(0.5, K + 0.5, by = 1))
  F_d     <- pgamma(pmax(edges_d, 0),
                    shape = par_d$shape, scale = par_d$scale)
  w_d     <- diff(F_d); w_d <- pmax(w_d, 0) / sum(pmax(w_d, 0))
  list(d = w_d, G = K, D = length(w_d))
}

# ============================================================================
# Data prep (idéntica a run_A_fsv.R)
# ============================================================================
wwtp_shp_batch <- sf::read_sf(shp_batch) |>
  sf::st_transform(32615) |>
  mutate(BatchName = toupper(trimws(BatchName)))
batch_key_shp <- wwtp_shp_batch |>
  sf::st_drop_geometry() |> distinct(BatchName) |> arrange(BatchName) |>
  mutate(Batch = paste0("WWTP", row_number()))
wwtp_shp_batch <- wwtp_shp_batch |> left_join(batch_key_shp, by = "BatchName")
wwtp_shp_focus <- wwtp_shp_batch |>
  filter(Batch %in% Batch_focus) |>
  arrange(match(Batch, Batch_focus))

ct_raw <- sf::st_drop_geometry(sf::read_sf(tracts_shp_path)) |>
  mutate(BatchName = toupper(trimws(BatchName)),
         MOE24_5YR = suppressWarnings(as.numeric(MOE24_5YR)))
ct_alloc <- ct_raw |> left_join(batch_key_shp, by = "BatchName") |>
  filter(!is.na(Batch)) |>
  mutate(pop_alloc = Percent_ce * POP24_5YR,
         moe_alloc = Percent_ce * MOE24_5YR)
ct_pop_by_batch <- ct_alloc |>
  group_by(Batch) |>
  summarise(N_acs = sum(pop_alloc, na.rm = TRUE),
            N_moe_acs = sqrt(sum(moe_alloc^2, na.rm = TRUE)),
            .groups = "drop")

batch_crosswalk <- tibble(
  batch_code = c("B69","BAS","BIA","BKW","BNO","BNE","BNW",
                 "BSB","BSO","BSE","BWE","BWS","BKB","BSW","BWB"),
  BatchName  = c("69TH STREET","ALMEDA SIMS","BATCH INTERCONTINENTAL",
                 "BATCH KINGWOOD","BATCH NORTH","BATCH NORTHEAST",
                 "BATCH NORTHWEST","BATCH SIMS BAYOU","BATCH SOUTH",
                 "BATCH SOUTHEAST","BATCH WEST","BATCH WEST SOUTH",
                 "KEEGANS BAYOU","SOUTHWEST","WILLOWBROOK")
)
ww_input <- read_csv(raw_csv_path, show_col_types = FALSE) |>
  mutate(WWTP = toupper(trimws(WWTP)),
         tn   = toupper(trimws(tn)),
         BATCH = toupper(trimws(BATCH))) |>
  filter(WWTP != "SS", tn == "RSV") |>
  transmute(batch_code = BATCH, WWTP, n_plants = Bnum,
            Date = as.Date(date), Pathogen = tn,
            copies_L, Flow = flow_interp) |>
  arrange(batch_code, Date) |>
  group_by(batch_code) |>
  mutate(median_flow = median(Flow, na.rm = TRUE)) |>
  ungroup() |>
  mutate(copies_day = (copies_L * median_flow) / 1e9) |>
  left_join(batch_crosswalk, by = "batch_code") |>
  left_join(batch_key_shp, by = "BatchName")
ww <- ww_input |> filter(Batch %in% Batch_focus) |> arrange(Batch, Date)

N_plantas <- ct_pop_by_batch |>
  filter(Batch %in% Batch_focus) |>
  arrange(match(Batch, Batch_focus)) |>
  { \(d) setNames(d$N_acs, d$Batch) }()
N_houston <- sum(ct_pop_by_batch$N_acs, na.rm = TRUE)
It_init_vec <- pmax(round(casos_houston_semana *
                          N_plantas[Batch_focus] / N_houston), 10)

ww_batch <- ww |>
  select(Batch, n_plants, Date, Pathogen, copies_day, flow_batch = Flow)
cutoff_date <- ww_batch |>
  summarise(d = dplyr::nth(sort(unique(Date)), burnin,
                           default = max(Date))) |>
  pull(d) |> max(cutoff_date_min)
ww_rsv_batch <- ww_batch |>
  filter(Date > cutoff_date) |>
  mutate(copies_day = ifelse(is.na(copies_day) | copies_day <= 0,
                              NA_real_, copies_day),
         Missing = is.na(copies_day))
make_wide_matrix <- function(df, value_col, plants_order) {
  df |> mutate(.val = {{ value_col }}) |>
    select(Date, WWTP = Batch, .val) |>
    dplyr::group_by(Date, WWTP) |>
    dplyr::summarise(.val = if (all(is.na(.val))) NA_real_ else mean(.val, na.rm = TRUE),
                     .groups = "drop") |>
    pivot_wider(names_from = WWTP, values_from = .val, values_fill = NA_real_) |>
    arrange(Date) |>
    { \(w) {
      missing_cols <- setdiff(plants_order, names(w))
      for (m in missing_cols) w[[m]] <- NA_real_
      w |> select(all_of(c("Date", plants_order)))
    } }()
}
y_wide_batch <- make_wide_matrix(ww_rsv_batch, copies_day, Batch_focus)
y_batch      <- as.matrix(y_wide_batch[, -1, drop = FALSE])
rownames(y_batch) <- as.character(y_wide_batch$Date)

P_focus <- ncol(y_batch); T_focus <- nrow(y_batch)
cat(sprintf("Dataset: %d weeks x %d plants\n", T_focus, P_focus))

# ============================================================================
# Kernels + W (sin MARSS)
# ============================================================================
ks <- make_shedding_kernel(mean_shed_days, sd_shed_days, K = 2)
d_vec <- as.numeric(ks$d); G <- ks$G; D <- ks$D
par_g <- to_gamma(prior_mean_gen_days, prior_sd_gen_days)
edges <- seq(1.5, G + 0.5, by = 1)
g_raw <- pgamma(edges, shape = par_g$shape, scale = par_g$scale) -
         pgamma(edges - 1, shape = par_g$shape, scale = par_g$scale)
g_raw <- pmax(g_raw, 0); g_vec <- g_raw / sum(g_raw)

# --- HAUSDORFF ABLATION: load precomputed W_HD.rds ------------------------
# The PEHD construction (region_distance_matrix(method="population", ...)
# followed by RBF kernel + row-normalize) has been replaced by loading the
# precomputed Hausdorff spatial weight matrix from disk. To regenerate it,
# run paper_pipeline/01_W_matrix/compute_W_matrices.R.
W_HD_path <- file.path(proj_root, "Modelos AR con rho est. y rho=0",
                         "paper_pipeline", "results", "W_HD.rds")
stopifnot(file.exists(W_HD_path))
W_HD <- readRDS(W_HD_path)
W_geo  <- W_HD[Batch_focus, Batch_focus, drop = FALSE]
cat(sprintf("Loaded W_HD from %s\n", W_HD_path))
cat(sprintf("W_HD dim: %d x %d, row-sums OK: %s\n",
            nrow(W_geo), ncol(W_geo),
            all(abs(rowSums(W_geo) - 1) < 1e-10)))

# ============================================================================
# Stan data — FSV con rho ESTIMADO (sin rho en data, sin sigma_obs / use_filter)
# ============================================================================
y_PT <- t(y_batch)
P <- length(Batch_focus); T <- ncol(y_PT)

y_obs_mat <- matrix(as.integer(is.finite(y_PT) & y_PT > 0),
                    nrow = nrow(y_PT), ncol = ncol(y_PT))
y_filled  <- y_PT
y_filled[y_obs_mat == 0] <- 1.0

cat(sprintf("Likelihood activa en %d / %d cells\n",
            sum(y_obs_mat), prod(dim(y_obs_mat))))

fits_dir <- file.path(proj_root, "outputs", "stan_fits_fsv_rhoest_HD")
dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
stan_path <- file.path(proj_root, "stan",
                       "renewal_A_fsv_est.stan")
stopifnot(file.exists(stan_path))

stan_data <- list(
  P = P, T = T, G = G, D = D,
  g = as.numeric(g_vec), d = as.numeric(d_vec),
  W = W_geo,
  # NO rho aqui — es parametro ahora
  It_init = as.numeric(It_init_vec),
  y = y_filled, y_obs = y_obs_mat,
  mean_gen_wk = prior_mean_gen_days / 7,
  sd_gen_wk   = prior_sd_gen_days   / 7,
  # NUEVO: hyperparametros del prior Beta(alpha_rho, beta_rho)
  alpha_rho = ALPHA_RHO,
  beta_rho  = BETA_RHO
)
init_list <- lapply(seq_len(n_chains), function(ch) list(
  log_beta = rep(0, P),
  sigma_obs_pl = rep(0.5, P),
  lam = matrix(0.1, P, T),
  It = matrix(rep(It_init_vec, T), P, T),
  rho = 0.5    # arrancar en el centro del prior
))

cat("\nCompiling Stan model...\n")
mod <- cmdstanr::cmdstan_model(stan_path, compile = TRUE,
                               cpp_options = list(stan_threads = TRUE))
cat("Compile done. Exe:", mod$exe_file(), "\n")

cat(sprintf("\n=== Sampling: %s ===\n", JOB_NAME))
t0 <- Sys.time()
fit <- mod$sample(
  data = stan_data, init = init_list,
  chains = n_chains, parallel_chains = n_chains,
  threads_per_chain = 1,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
  refresh = 500,
  seed = 1729    # seed fijo (no depende de rho)
)
cat(sprintf("\nSample done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("Extracting draws (incluye rho)...\n")
draws_arr <- fit$draws(
  variables = c("log_beta", "beta", "sigma_obs_pl", "Rt", "It",
                "rho", "log_lik", "lp__"),
  format = "draws_array"
)
diag_sum  <- fit$diagnostic_summary()
summary_tbl <- fit$summary(
  variables = c("log_beta", "beta", "sigma_obs_pl", "Rt", "It", "rho", "lp__"),
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail,
  mean, sd,
  q5 = ~ quantile(.x, 0.05),
  q50 = ~ quantile(.x, 0.50),
  q95 = ~ quantile(.x, 0.95)
)

# Print rho posterior summary front-and-center
rho_sm <- summary_tbl[summary_tbl$variable == "rho", ]
cat("\n=== POSTERIOR DE RHO ===\n")
print(rho_sm)

save_object <- list(
  rho_prior = list(family = "beta", alpha = ALPHA_RHO, beta = BETA_RHO),
  rho_posterior_summary = rho_sm,
  model_variant = "freesigma_vec_rhoest",
  draws = draws_arr,
  diagnostic_summary = diag_sum,
  summary = summary_tbl,
  y_raw = y_PT,
  y_obs_mask = y_obs_mat,
  d_vec = d_vec,
  meta = list(
    n_chains = n_chains, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    log_beta_prior  = "N(0, 4)",
    sigma_obs_prior = "Half-Normal(0, 2) por planta",
    rho_prior       = sprintf("Beta(%g, %g)", ALPHA_RHO, BETA_RHO),
    save_time = Sys.time()
  )
)

final_path <- file.path(fits_dir, sprintf("%s.rds", JOB_NAME))
tmp_path <- file.path(Sys.getenv("TMPDIR", tempdir()),
                      sprintf("%s_%d.rds", JOB_NAME, Sys.getpid()))
saveRDS(save_object, tmp_path, compress = "xz")
test_read <- tryCatch(readRDS(tmp_path), error = function(e) NULL)
if (is.null(test_read)) stop("Tmp corrupto")
cat(sprintf("Tmp OK (%.1f MB)\n", file.size(tmp_path) / 1e6))
ok <- file.copy(tmp_path, final_path, overwrite = TRUE)
if (!ok) stop("file.copy fallo")
file.remove(tmp_path)
final_read <- tryCatch(readRDS(final_path), error = function(e) NULL)
if (is.null(final_read)) stop("Final corrupto")
cat(sprintf("Final OK (%.1f MB)\n", file.size(final_path) / 1e6))

cat(sprintf("\n=== %s DONE ===\n", JOB_NAME))
