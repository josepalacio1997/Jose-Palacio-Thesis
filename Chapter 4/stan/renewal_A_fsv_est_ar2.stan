// ============================================================================
// rsv_renewal_freesigma_vec_rhoest_ar2.stan
//
// AR(2) variant of rsv_renewal_freesigma_vec_rhoest.stan.
// Cambio UNICO vs el rho-est AR(1): la mean del growth rate lambda
// pasa de depender de UN lag a depender de DOS lags con coeficientes
// FIJOS (phi1, phi2) = (2, -1) — los mismos que MARSS usa sobre el nivel
// log10(RNA) en Ensor et al. 2025.
//
//   AR(1) (rho-est original):
//     mu_lam[p,t] = (1-rho) * lam[p,t-1] + rho * sum_j W[p,j] * lam[j,t-1]
//
//   AR(2) (esta version):
//     spatial_t1 = (1-rho)*lam[p,t-1] + rho*sum_j W[p,j]*lam[j,t-1]
//     spatial_t2 = (1-rho)*lam[p,t-2] + rho*sum_j W[p,j]*lam[j,t-2]
//     mu_lam[p,t] = 2 * spatial_t1  -  1 * spatial_t2
//
// Justificacion: MARSS hace l_t = 2*l_{t-1} - l_{t-2} + w_t, es decir
// AR(2) en el nivel con coefs (2, -1). Aplicar la misma forma a lambda
// aqui pone el modelo en correspondencia estructural con MARSS de manera
// directa, sin agregar parametros.
//
// Initial conditions: lambda[p,1] y lambda[p,2] ambos ~ N(0, 0.5)
// (mismo prior que el AR(1) usaba para lambda[p,1]).
//
// Todo el resto del modelo (priors, sigma_obs_pl vector[P], beta, It,
// shedding kernel, lognormal observation) queda IDENTICO al rho-est AR(1).
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
  vector<lower=1>[P] It_init;
  matrix<lower=0>[P, T] y;
  array[P, T] int<lower=0, upper=1> y_obs;
  real<lower=0> mean_gen_wk;
  real<lower=0> sd_gen_wk;
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
}

parameters {
  vector<lower=-20, upper=20>[P] log_beta;
  matrix<lower=-20, upper=20>[P, T] lam;
  matrix<lower=1e-6, upper=1e10>[P, T] It;
  vector<lower=0, upper=10>[P] sigma_obs_pl;
  real<lower=0, upper=1> rho;
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
  // ---- PRIORS (unchanged) ---------------------------------------------------
  log_beta ~ normal(0, 4);
  sigma_obs_pl ~ normal(0, 2);
  rho ~ beta(alpha_rho, beta_rho);

  // ---- AR(2) on lambda with fixed (phi1, phi2) = (2, -1) -------------------
  // Two anchors instead of one (mismo prior N(0, 0.5) que el AR(1) usaba)
  for (p in 1:P) lam[p, 1] ~ normal(0, 0.5);
  if (T >= 2) for (p in 1:P) lam[p, 2] ~ normal(0, 0.5);

  // Recursion starts at t = 3 because we need lam[, t-1] AND lam[, t-2]
  for (t in 3:T) {
    for (p in 1:P) {
      // Spatial mix at t-1
      real nbr_t1 = 0;
      for (j in 1:P) nbr_t1 += W[p, j] * lam[j, t - 1];
      real spatial_t1 = (1 - rho) * lam[p, t - 1] + rho * nbr_t1;

      // Spatial mix at t-2
      real nbr_t2 = 0;
      for (j in 1:P) nbr_t2 += W[p, j] * lam[j, t - 2];
      real spatial_t2 = (1 - rho) * lam[p, t - 2] + rho * nbr_t2;

      // AR(2) with fixed phi1 = 2, phi2 = -1 (MARSS-equivalent on level)
      real mu_lam = 2.0 * spatial_t1 - 1.0 * spatial_t2;
      lam[p, t] ~ normal(mu_lam, 0.5);
    }
  }

  // ---- LATENT It (Poisson -> LogNormal moment-matched relaxation) ----------
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

  // ---- OBSERVATION (lognormal, sigma_obs_pl free per plant) ----------------
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(sigma_obs_pl[p]);
        y[p, t] ~ lognormal(mu_obs, sigma_obs_pl[p]);
      }
    }
  }
}

generated quantities {
  matrix[P, T] log_lik;
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(sigma_obs_pl[p]);
        log_lik[p, t] = lognormal_lpdf(y[p, t] | mu_obs, sigma_obs_pl[p]);
      } else {
        log_lik[p, t] = 0;
      }
    }
  }
}
