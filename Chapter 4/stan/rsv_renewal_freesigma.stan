// ============================================================================
// rsv_renewal_freesigma.stan
// Variante del modelo principal SIN MARSS Stage 1:
//   - sigma_obs NO se pasa como dato; se estima como UN solo parametro
//     escalar compartido entre todas las plantas y todas las semanas.
//   - log_beta ~ Normal(0, 4): prior vago (no informativo en el rango
//     biologicamente plausible).
//   - Half-Normal(0, 2) sobre sigma_obs_pl: prior weakly-informative
//     estandar para variance components (Gelman 2006, Bayesian Analysis).
// ============================================================================

functions {
  real euler_lotka_rt(real lambda, real mean_gen, real sd_gen) {
    real arg = 1 + (square(sd_gen) / mean_gen) * lambda;
    if (arg <= 0) return 1e-6;
    real exponent = square(mean_gen) / square(sd_gen);
    return pow(arg, exponent);
  }
}

data {
  int<lower=1> P;
  int<lower=1> T;
  int<lower=1> G;
  int<lower=1> D;
  vector<lower=0>[G] g;
  vector<lower=0>[D] d;
  matrix[P, P] W;
  real<lower=0, upper=1> rho;
  vector<lower=1>[P] It_init;
  matrix<lower=0>[P, T] y;
  array[P, T] int<lower=0, upper=1> y_obs;
  real<lower=0> mean_gen_wk;
  real<lower=0> sd_gen_wk;
}

parameters {
  vector<lower=-20, upper=20>[P] log_beta;
  matrix<lower=-20, upper=20>[P, T] lam;
  matrix<lower=1e-6, upper=1e10>[P, T] It;

  // SINGLE shared observation noise (scalar, no indexing).
  // Upper bound 10 keeps sampler from drifting; in practice sigma converges
  // to ~0.6-1.5 based on noise-floor analysis of the data.
  real<lower=0, upper=10> sigma_obs_pl;
}

transformed parameters {
  vector<lower=0>[P] beta = exp(log_beta);
  matrix[P, T] Rt;
  matrix<lower=0>[P, T] iota;
  matrix<lower=0>[P, T] I_mix;
  matrix<lower=0>[P, T] psi;

  for (p in 1:P) {
    Rt[p, 1] = euler_lotka_rt(lam[p, 1], mean_gen_wk, sd_gen_wk);
    iota[p, 1] = It_init[p];
  }
  for (p in 1:P) {
    real sp = 0;
    for (j in 1:P) sp += W[p, j] * It[j, 1];
    I_mix[p, 1] = (1 - rho) * It[p, 1] + rho * sp;
  }

  for (t in 2:T) {
    for (p in 1:P) Rt[p, t] = euler_lotka_rt(lam[p, t], mean_gen_wk, sd_gen_wk);
    for (p in 1:P) {
      real conv = 0;
      for (k in 1:G) {
        int s = t - k;
        if (s >= 1) conv += g[k] * I_mix[p, s];
      }
      iota[p, t] = Rt[p, t] * fmax(conv, 1e-10);
    }
    for (p in 1:P) {
      real sp = 0;
      for (j in 1:P) sp += W[p, j] * It[j, t];
      I_mix[p, t] = (1 - rho) * It[p, t] + rho * sp;
    }
  }

  for (p in 1:P) {
    for (t in 1:T) {
      real sh = 0;
      for (k in 0:(D - 1)) {
        int s = t - k;
        if (s >= 1) sh += d[k + 1] * It[p, s];
      }
      psi[p, t] = beta[p] * fmax(sh, 1e-10);
    }
  }
}

model {
  // PRIORS ------------------------------------------------------------------
  // Weakly-informative on log_beta: ~5 orders of magnitude on beta,
  // essentially flat over the biologically plausible range.
  log_beta ~ normal(0, 4);

  // Half-Normal(0, 2) on the shared observation noise. Lower=0 + this
  // prior puts mass in [0, ~4] with mode at 0, decreasing tail.
  // Reference: Gelman (2006) "Prior distributions for variance parameters
  // in hierarchical models." Bayesian Analysis 1(3): 515-533.
  sigma_obs_pl ~ normal(0, 2);

  // Latent log-growth-rate AR with spatial coupling
  for (p in 1:P) lam[p, 1] ~ normal(0, 0.5);
  for (t in 2:T) {
    for (p in 1:P) {
      real mu_lam = 0;
      for (j in 1:P) mu_lam += W[p, j] * lam[j, t - 1];
      mu_lam = (1 - rho) * lam[p, t - 1] + rho * mu_lam;
      lam[p, t] ~ normal(mu_lam, 0.5);
    }
  }

  // Poisson -> LogNormal moment-matched relaxation, clamped at both ends.
  // Floor at 10 prevents log(lam) - 0.5/lam from -inf when iota -> 0.
  // Cap at 1e8 prevents +inf during wild NUTS warmup proposals.
  for (p in 1:P) {
    real lam_p = fmin(fmax(It_init[p], 10.0), 1e8);
    real mu_ln = log(lam_p) - 0.5 / lam_p;
    real sd_ln = inv_sqrt(lam_p);
    target += lognormal_lpdf(It[p, 1] | mu_ln, sd_ln);
  }
  for (t in 2:T) {
    for (p in 1:P) {
      real lam_p = fmin(fmax(iota[p, t], 10.0), 1e8);
      real mu_ln = log(lam_p) - 0.5 / lam_p;
      real sd_ln = inv_sqrt(lam_p);
      target += lognormal_lpdf(It[p, t] | mu_ln, sd_ln);
    }
  }

  // Observation model: shared sigma_obs_pl (scalar) across all (p, t).
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(sigma_obs_pl);
        y[p, t] ~ lognormal(mu_obs, sigma_obs_pl);
      }
    }
  }
}

generated quantities {
  matrix[P, T] log_lik;
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(sigma_obs_pl);
        log_lik[p, t] = lognormal_lpdf(y[p, t] | mu_obs, sigma_obs_pl);
      } else {
        log_lik[p, t] = 0;
      }
    }
  }
}
