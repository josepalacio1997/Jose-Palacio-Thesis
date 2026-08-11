# Extract rho posterior summaries for Filter and Data variants
suppressPackageStartupMessages({
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

BASE <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0/fits"

fits <- list(
  Filter = readRDS(file.path(BASE, "A_ssm_est.rds")),
  Data   = readRDS(file.path(BASE, "A_direct_est.rds"))
)

summarize_rho <- function(obj, label) {
  # Cross-check: use pre-computed summary (fast) AND recompute q025/q975 from draws
  s <- obj$rho_posterior_summary
  d <- as_draws_matrix(obj$draws)
  r <- as.numeric(d[, "rho"])
  data.frame(
    variant = label,
    mean    = s$mean,
    sd      = s$sd,
    q05     = s$`5%`,
    q95     = s$`95%`,
    q025    = quantile(r, 0.025),
    q975    = quantile(r, 0.975),
    rhat    = s$rhat,
    ess_bulk = s$ess_bulk,
    row.names = NULL
  )
}

result <- do.call(rbind, Map(summarize_rho, fits, names(fits)))
print(result, row.names = FALSE, digits = 4)

out_dir <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0/results_modelA_AR1"
write.csv(result, file.path(out_dir, "rho_summary.csv"), row.names = FALSE)
cat("\n  -> rho_summary.csv written to results_modelA_AR1/\n")

# ---- Posterior density plot with 95% CrI shaded ----
draws_long <- do.call(rbind, lapply(names(fits), function(nm) {
  d <- as_draws_matrix(fits[[nm]]$draws)
  data.frame(variant = nm, rho = as.numeric(d[, "rho"]))
}))

crI_df <- do.call(rbind, lapply(names(fits), function(nm) {
  d <- as_draws_matrix(fits[[nm]]$draws)
  r <- as.numeric(d[, "rho"])
  data.frame(variant = nm,
             q025 = quantile(r, 0.025),
             q975 = quantile(r, 0.975),
             mean = mean(r))
}))

p_rho <- ggplot(draws_long, aes(x = rho, fill = variant, color = variant)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(data = crI_df, aes(xintercept = mean, color = variant),
             linetype = "dashed", linewidth = 0.5) +
  geom_segment(data = crI_df,
               aes(x = q025, xend = q975, y = 0, yend = 0, color = variant),
               linewidth = 2, alpha = 0.8) +
  scale_fill_manual(values  = c(Filter = "#2C7BB6", Data = "#E67E22")) +
  scale_color_manual(values = c(Filter = "#2C7BB6", Data = "#E67E22")) +
  labs(x = expression(rho),
       y = "Posterior density",
       title = expression("Posterior of "*rho*" with 95% credible intervals"),
       fill = "Variant", color = "Variant") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")

ggsave(file.path(out_dir, "rho_posterior_95CrI.pdf"),
       p_rho, width = 8, height = 5)
cat("  -> rho_posterior_95CrI.pdf written to results_modelA_AR1/\n")
