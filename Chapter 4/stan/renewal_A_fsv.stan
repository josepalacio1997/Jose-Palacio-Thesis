// ============================================================================
// rsv_renewal_freesigma_vec.stan
// Variante de rsv_renewal_freesigma.stan con sigma_obs_pl POR PLANTA (vector[P]).
//
// Cambios vs rsv_renewal_freesigma.stan:
//   - sigma_obs_pl: real -> vector[P]
//   - Half-Normal(0, 2) prior aplicado componentwise (Stan auto-vectoriza)
//   - mu_obs y likelihood usan sigma_obs_pl[p] en lugar del escalar
//
// Motivacion: el modelo escalar produce coverage 95% empirico = 50% por
// heterogeneidad real entre plantas (noise floor va de 0.47 a 0.98).
// Vectorizar permite a cada planta tener su propio sigma_p y reduce
// la sobreconfianza en el observation model.
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

  // CAMBIO PRINCIPAL: vector[P] en vez de real
  // Cada planta tiene su propio sigma_obs. Half-Normal(0, 2) prior.
  vector<lower=0, upper=10>[P] sigma_obs_pl;
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
  log_beta ~ normal(0, 4);

  // Half-Normal(0, 2) componentwise sobre el vector[P] sigma_obs_pl.
  // Stan auto-vectoriza: equivale a sum_{p=1}^P normal_lpdf(sigma_obs_pl[p] | 0, 2).
  sigma_obs_pl ~ normal(0, 2);

  for (p in 1:P) lam[p, 1] ~ normal(0, 0.5);
  for (t in 2:T) {
    for (p in 1:P) {
      real mu_lam = 0;
      for (j in 1:P) mu_lam += W[p, j] * lam[j, t - 1];
      mu_lam = (1 - rho) * lam[p, t - 1] + rho * mu_lam;
      lam[p, t] ~ normal(mu_lam, 0.5);
    }
  }

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

  // Observation: sigma_obs_pl[p] (por planta) en vez del escalar.
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
