# ============================================================================
# run_A_marss_est.R
#
# In-sample fit del modelo MARSS-hybrid (Ensor-style) con RHO ESTIMADO.
# Una sola corrida (sin loop rho), prior Beta(alpha_rho, beta_rho).
# Llamada: Rscript run_A_marss_est.R [ALPHA_RHO] [BETA_RHO] [WARMUP] [SAMPLING]
#   ej. Rscript run_A_marss_est.R                  # defaults: Beta(2,2)
#   ej. Rscript run_A_marss_est.R 2 2 4000 6000    # explicit
#   ej. Rscript run_A_marss_est.R 1 1              # Uniform(0,1)
#
# Stan: stan/renewal_A_marss_est.stan (use_filter=1, sigma_obs=MARSS Stage 1)
# Output: outputs/stan_fits_marss_rhoest/rho_estimated.rds
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
ALPHA_RHO    <- if (length(args) >= 1) as.numeric(args[1]) else 2.0
BETA_RHO     <- if (length(args) >= 2) as.numeric(args[2]) else 2.0
arg_warmup   <- if (length(args) >= 3) suppressWarnings(as.integer(args[3])) else NA_integer_
arg_sampling <- if (length(args) >= 4) suppressWarnings(as.integer(args[4])) else NA_integer_
stopifnot(ALPHA_RHO > 0, BETA_RHO > 0)
JOB_NAME <- "rho_estimated"
cat(sprintf("\n=== IN-SAMPLE MARSS (rho ESTIMATED, Beta(%g, %g)) ===\n",
            ALPHA_RHO, BETA_RHO))

suppressPackageStartupMessages({
  library(tidyverse); library(lubridate); library(sf); library(MARSS)
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
# Data prep (idéntica a run_A_marss.R)
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
# MARSS Stage 1 — produce y_filt y sigma_filt para CADA (p, t)
# ============================================================================
ln10 <- log(10)
Y_log10_TP <- log10(y_batch)
Y_log10_TP[!is.finite(Y_log10_TP)] <- NA

ssm_spec_one_plant <- list(
  B  = matrix(c(2, 1, -1, 0), 2, 2, byrow = TRUE),
  U  = "zero", Q = "diagonal and equal", V0 = "identity", x0 = "equal",
  R  = "diagonal and equal", A = "zero", Z = matrix(c(1, 0), 1, 2)
)
y_filt     <- matrix(NA_real_, P_focus, T_focus,
                     dimnames = list(Batch_focus, NULL))
sigma_filt <- matrix(NA_real_, P_focus, T_focus,
                     dimnames = list(Batch_focus, NULL))
cat("Running MARSS Stage 1...\n")
for (p in seq_len(P_focus)) {
  yp_log10 <- Y_log10_TP[, p]
  if (sum(!is.na(yp_log10)) < 5) {
    y_filt[p, ] <- 1; sigma_filt[p, ] <- 1; next
  }
  fit_p <- MARSS::MARSS(matrix(yp_log10, nrow = 1),
                        model = ssm_spec_one_plant,
                        method = "BFGS", silent = TRUE,
                        control = list(maxit = 200))
  kfas_p <- MARSS::MARSSkfas(fit_p)
  xtt_p <- kfas_p$xtt; Vtt_p <- kfas_p$Vtt
  Mp <- coef(fit_p, type = "matrix")
  Zp <- as.matrix(Mp$Z); Rp <- as.matrix(Mp$R)
  mu_log10_t <- as.numeric(Zp %*% xtt_p)
  V_obs_t <- vapply(seq_len(T_focus),
                    function(tt) as.numeric(Zp %*% Vtt_p[, , tt] %*% t(Zp)),
                    numeric(1))
  V_obs_t <- pmax(V_obs_t, 0)
  mu_ln_t  <- ln10 * mu_log10_t
  var_ln_t <- (ln10^2) * V_obs_t
  y_filt[p, ] <- exp(mu_ln_t + 0.5 * var_ln_t)
  sigma_filt[p, ] <- sqrt(var_ln_t)
}
sigma_filt <- pmin(pmax(sigma_filt, 1e-6), 3.0)
cat("MARSS done. y_filt range: [", signif(min(y_filt), 3), ",",
    signif(max(y_filt), 3), "]\n", sep = "")

# ============================================================================
# Kernels + W
# ============================================================================
ks <- make_shedding_kernel(mean_shed_days, sd_shed_days, K = 2)
d_vec <- as.numeric(ks$d); G <- ks$G; D <- ks$D
par_g <- to_gamma(prior_mean_gen_days, prior_sd_gen_days)
edges <- seq(1.5, G + 0.5, by = 1)
g_raw <- pgamma(edges, shape = par_g$shape, scale = par_g$scale) -
         pgamma(edges - 1, shape = par_g$shape, scale = par_g$scale)
g_raw <- pmax(g_raw, 0); g_vec <- g_raw / sum(g_raw)

tracts_raw <- sf::st_read(tracts_shp_path, quiet = TRUE) |>
  sf::st_make_valid() |>
  dplyr::filter(sf::st_is_valid(geometry), !sf::st_is_empty(geometry))
tracts_valid <- tracts_raw |>
  dplyr::distinct(GEOID, .keep_all = TRUE) |>
  sf::st_transform(32615)
d_pop_50 <- region_distance_matrix(
  regions = wwtp_shp_focus, label_col = "BatchName",
  method = "population", tracts = tracts_valid, pop_col = "POP24_5YR",
  target_fraction = 0.5, people_per_dot = 200, units = "km")
make_W_from_distance <- function(d, sigma = NULL) {
  d_dim <- dim(d); d_names <- dimnames(d)
  d <- structure(as.numeric(unclass(d)), dim = d_dim, dimnames = d_names)
  off <- d; diag(off) <- NA
  if (is.null(sigma)) sigma <- sqrt(median(off^2, na.rm = TRUE) / 2)
  K <- exp(-d^2 / (2 * sigma^2)); diag(K) <- 0
  W <- K / rowSums(K); W
}
W_pop50 <- make_W_from_distance(d_pop_50)
bn_to_wwtp <- setNames(batch_key_shp$Batch, batch_key_shp$BatchName)
rownames(W_pop50) <- bn_to_wwtp[rownames(W_pop50)]
colnames(W_pop50) <- bn_to_wwtp[colnames(W_pop50)]
W_geo <- W_pop50[Batch_focus, Batch_focus, drop = FALSE]

# ============================================================================
# Stan data — MARSS rho ESTIMADO
# ============================================================================
y_PT <- t(y_batch)
P <- length(Batch_focus); T <- ncol(y_PT)

y_filled  <- y_filt                                          # MARSS-imputado
y_obs_mat <- matrix(1L, nrow = nrow(y_PT), ncol = ncol(y_PT))  # all-ones
y_obs_raw <- matrix(as.integer(is.finite(y_PT) & y_PT > 0),
                    nrow = nrow(y_PT), ncol = ncol(y_PT))
cat(sprintf("Likelihood activa en %d / %d cells (Ensor-style)\n",
            sum(y_obs_mat), prod(dim(y_obs_mat))))
cat(sprintf("  Detalle: %d cells con y_raw real, %d cells MARSS-imputado\n",
            sum(y_obs_raw), sum(y_obs_mat) - sum(y_obs_raw)))

fits_dir <- file.path(proj_root, "outputs", "stan_fits_marss_rhoest")
dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
stan_path <- file.path(proj_root, "stan", "renewal_A_marss_est.stan")
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
  use_filter = 1L,
  sigma_obs = sigma_filt,
  # NUEVO: hyperparametros del prior Beta(alpha_rho, beta_rho)
  alpha_rho = ALPHA_RHO,
  beta_rho  = BETA_RHO
)
init_list <- lapply(seq_len(n_chains), function(ch) list(
  log_beta = rep(0, P), sigma_obs_pl = rep(0.5, P),
  lam = matrix(0.1, P, T), It = matrix(rep(It_init_vec, T), P, T),
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

rho_sm <- summary_tbl[summary_tbl$variable == "rho", ]
cat("\n=== POSTERIOR DE RHO ===\n")
print(rho_sm)

save_object <- list(
  rho_prior = list(family = "beta", alpha = ALPHA_RHO, beta = BETA_RHO),
  rho_posterior_summary = rho_sm,
  model_variant = "MARSS_hybrid_rhoest",
  draws = draws_arr,
  diagnostic_summary = diag_sum,
  summary = summary_tbl,
  y_raw = y_PT,
  y_filt = y_filt,
  sigma_filt = sigma_filt,
  y_obs_mask = y_obs_raw,
  d_vec = d_vec,
  meta = list(
    n_chains = n_chains, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    log_beta_prior  = "N(0, 4)",
    sigma_obs_prior = "lognormal(log(0.5), sqrt(10))  (no usado: use_filter=1)",
    sigma_obs_data  = "MARSS Stage 1 per-(p,t) variance",
    rho_prior       = sprintf("Beta(%g, %g)", ALPHA_RHO, BETA_RHO),
    save_time = Sys.time()
  )
)

final_path <- file.path(fits_dir, sprintf("%s.rds", JOB_NAME))
# ---- LUSTRE-SAFE SAVE v3: direct saveRDS + sync + verify ----
# tmp + file.copy genera .rds truncados en Lustre porque file.copy returna
# TRUE sin que el write async haya commiteado a los OSTs. v3 elimina el
# tmp y va directo a Lustre con gc() previo + system("sync") + verify.
cat(sprintf("Save object size in R memory: ~%.1f MB\n",
            as.numeric(object.size(save_object)) / 1e6))

cat("Pre-save: rm(fit) + gc() para liberar memoria...\n")
rm(fit)
invisible(gc(verbose = FALSE))

cat(sprintf("Step 1/4: saveRDS directo a Lustre (%s) ...\n", final_path))
saveRDS(save_object, final_path, compress = "gzip")
cat("Step 1/4 OK: saveRDS returned\n")

cat("Step 2/4: shell sync para forzar fsync a Lustre OSTs ...\n")
sync_ret <- system("sync")
cat(sprintf("Step 2/4 OK: sync returned %d\n", sync_ret))

cat("Step 3/4: Sys.sleep(30) buffer + verificar tamano ...\n")
Sys.sleep(30)
sz_mb <- file.size(final_path) / 1e6
cat(sprintf("Step 3/4: file size en Lustre = %.2f MB\n", sz_mb))
if (sz_mb < 50) stop(sprintf("Step 3/4 FAIL: file too small (%.2f MB)", sz_mb))

cat("Step 4/4: readRDS verify (con retry de 60s) ...\n")
final_read <- tryCatch(readRDS(final_path),
                       error = function(e) {
                         cat("First readRDS failed:", conditionMessage(e), "\n")
                         cat("Waiting 60s and retrying...\n")
                         Sys.sleep(60)
                         tryCatch(readRDS(final_path),
                                  error = function(e2) {
                                    cat("Second readRDS failed:",
                                        conditionMessage(e2), "\n")
                                    NULL
                                  })
                       })
if (is.null(final_read)) stop("Step 4/4 FAIL: file unreadable after retry")
cat(sprintf("ALL OK: file %.2f MB en Lustre, leible.\n", sz_mb))

cat(sprintf("\n=== %s DONE ===\n", JOB_NAME))
