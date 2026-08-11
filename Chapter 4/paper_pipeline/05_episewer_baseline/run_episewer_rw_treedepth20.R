# Variante de run_episewer_rw.R para plant_3 con max_treedepth=20
# (plant_3 saturó max_treedepth=15 default)

args <- commandArgs(trailingOnly = TRUE)
PLANT_IDX     <- as.integer(args[1])
iter_warmup   <- if (length(args) >= 2) as.integer(args[2]) else 5000
iter_sampling <- if (length(args) >= 3) as.integer(args[3]) else 5000

stopifnot(PLANT_IDX >= 1, PLANT_IDX <= 15)
plant_name <- sprintf("WWTP%d", PLANT_IDX)
cat(sprintf("\n=== EpiSewer RW (max_treedepth=20) for %s ===\n", plant_name))

suppressPackageStartupMessages({
  library(EpiSewer)
  library(dplyr); library(readr); library(lubridate); library(sf); library(cmdstanr)
  library(data.table)
})

set.seed(42)
proj_root <- normalizePath(".", mustWork = TRUE)

# ---- Data prep (FULL median, igual que run_episewer_rw.R) ------------------
raw_csv  <- file.path(proj_root, "data",
                      "all_ts_observed_ILI_20250526_pool.csv")
shp_path <- file.path(proj_root, "data", "WWTP_Batch",
                      "WWTP_061621_mergedWWTPs to Batch_updatedon07192024.shp")

shp <- sf::read_sf(shp_path) |> sf::st_drop_geometry() |>
  mutate(BatchName = toupper(trimws(BatchName))) |>
  distinct(BatchName) |> arrange(BatchName) |>
  mutate(Batch = paste0("WWTP", row_number()))
batch_name <- shp$BatchName[shp$Batch == plant_name]

batch_crosswalk <- tibble::tibble(
  batch_code = c("B69","BAS","BIA","BKW","BNO","BNE","BNW",
                 "BSB","BSO","BSE","BWE","BWS","BKB","BSW","BWB"),
  BatchName  = c("69TH STREET","ALMEDA SIMS","BATCH INTERCONTINENTAL",
                 "BATCH KINGWOOD","BATCH NORTH","BATCH NORTHEAST",
                 "BATCH NORTHWEST","BATCH SIMS BAYOU","BATCH SOUTH",
                 "BATCH SOUTHEAST","BATCH WEST","BATCH WEST SOUTH",
                 "KEEGANS BAYOU","SOUTHWEST","WILLOWBROOK")
)
batch_code <- batch_crosswalk$batch_code[batch_crosswalk$BatchName == batch_name]
cat(sprintf("BatchName: %s | batch_code: %s\n", batch_name, batch_code))

ww_full <- read_csv(raw_csv, show_col_types = FALSE) |>
  mutate(WWTP = toupper(trimws(WWTP)),
         tn   = toupper(trimws(tn)),
         BATCH= toupper(trimws(BATCH))) |>
  filter(WWTP != "SS", tn == "RSV", BATCH == batch_code) |>
  arrange(date)
median_flow <- median(ww_full$flow_interp, na.rm = TRUE)
cat(sprintf("median_flow FULL: %.2e L/day\n", median_flow))

ww_raw <- ww_full |> filter(date > as.Date("2023-01-02"))
ww <- ww_raw |> filter(!is.na(copies_L), copies_L > 0)
n_dates <- length(unique(ww$date))
if (n_dates < nrow(ww)) {
  ww <- ww |> group_by(date) |>
    summarise(copies_L = mean(copies_L, na.rm = TRUE),
              flow_interp = mean(flow_interp, na.rm = TRUE),
              .groups = "drop")
}

ww <- ww |> arrange(date)
ref_date <- as.Date("2020-01-06")
pseudo_dates <- ref_date + seq_len(nrow(ww)) - 1
ww_dt <- data.table::data.table(
  date = pseudo_dates,
  real_date = as.Date(ww$date),
  concentration = (ww$copies_L * median_flow) / 1e9,
  flow = 1,
  cases = NA_real_
)
cat(sprintf("n_obs: %d, y_raw range [%.2f, %.2f]\n",
            nrow(ww_dt), min(ww_dt$concentration), max(ww_dt$concentration)))

ww_sewer <- list(
  measurements = ww_dt[, .(date, concentration)],
  flows = ww_dt[, .(date, flow)],
  cases = ww_dt[, .(date, cases)]
)

# ---- Module configs (igual que original) -----------------------------------
ep_measurements <- EpiSewer::model_measurements(
  concentrations = EpiSewer::concentrations_observe(measurements = ww_sewer$measurements),
  noise = EpiSewer::noise_estimate(),
  LOD = EpiSewer::LOD_none()
)
ep_sampling <- EpiSewer::model_sampling(
  sample_effects = EpiSewer::sample_effects_none()
)
ep_sewage <- EpiSewer::model_sewage(
  flows = EpiSewer::flows_observe(flows = ww_sewer$flows),
  residence_dist = EpiSewer::residence_dist_assume(residence_dist = c(1))
)
ep_shedding <- EpiSewer::model_shedding(
  shedding_dist = EpiSewer::shedding_dist_assume(
    EpiSewer::get_discrete_gamma(gamma_mean = 4.6, gamma_sd = 2.0, maxX = 14),
    shedding_reference = "symptom_onset"),
  incubation_dist = EpiSewer::incubation_dist_assume(
    EpiSewer::get_discrete_gamma(gamma_mean = 5.0, gamma_sd = 2.0, maxX = 14)),
  load_per_case = EpiSewer::load_per_case_calibrate(cases = NULL, min_cases = 10),
  load_variation = EpiSewer::load_variation_estimate()
)
ep_infections <- EpiSewer::model_infections(
  generation_dist = EpiSewer::generation_dist_assume(
    EpiSewer::get_discrete_gamma_shifted(gamma_mean = 7.5, gamma_sd = 2.1, maxX = 21)),
  R = EpiSewer::R_estimate_rw(),
  seeding = EpiSewer::seeding_estimate_rw(),
  infection_noise = EpiSewer::infection_noise_estimate()
)
ep_forecast <- EpiSewer::model_forecast(
  horizon = EpiSewer::horizon_assume(horizon = 7)
)

# *** KEY CHANGE: max_treedepth = 20 (no 15) ***
ep_fit_opts <- EpiSewer::set_fit_opts(
  model = EpiSewer::model_stan_opts(package = "EpiSewer"),
  sampler = EpiSewer::sampler_stan_mcmc(
    chains = 4, parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = 0.99,
    max_treedepth = 20,             # <-- bump 15 -> 20
    seed = 1729
  )
)

ep_results_opts <- EpiSewer::set_results_opts(
  fitted = TRUE,
  summary_intervals = c(0.5, 0.95),
  samples_ndraws = 100
)

cat(sprintf("\nSampling %s with %d+%d iters, max_treedepth=20\n",
            plant_name, iter_warmup, iter_sampling))
t0 <- Sys.time()
ep_result <- EpiSewer::EpiSewer(
  measurements = ep_measurements,
  sampling     = ep_sampling,
  sewage       = ep_sewage,
  shedding     = ep_shedding,
  infections   = ep_infections,
  forecast     = ep_forecast,
  fit_opts     = ep_fit_opts,
  results_opts = ep_results_opts
)
cat(sprintf("\nEpiSewer done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

out_dir <- file.path(proj_root, "Modelos AR con rho est. y rho=0",
                     "episewer_rw_fits")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rds_path <- file.path(out_dir, sprintf("plant_%02d.rds", PLANT_IDX))

save_object <- list(
  plant_idx = PLANT_IDX, plant_name = plant_name, batch_name = batch_name,
  batch_code = batch_code, result = ep_result,
  input_data = data.frame(
    real_date = as.Date(ww$date),
    pseudo_date = pseudo_dates,
    copies_L = ww$copies_L,
    measurement = (ww$copies_L * median_flow) / 1e9
  ),
  meta = list(
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    units = sprintf("daily RW, %d+%d iters, max_treedepth=20", iter_warmup, iter_sampling),
    median_flow_constant = median_flow,
    save_time = Sys.time()
  )
)
saveRDS(save_object, rds_path, compress = "gzip")
cat(sprintf("Saved -> %s  (%.1f MB)\n",
            rds_path, file.size(rds_path) / 1e6))
cat("=== DONE ===\n")
