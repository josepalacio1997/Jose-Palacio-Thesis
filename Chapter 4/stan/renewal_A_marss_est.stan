// ============================================================================
// rsv_renewal_rhoest.stan
//
// Variante de rsv_renewal.stan (MARSS) donde rho se ESTIMA como parametro
// (en vez de pasarse como data y hacer sweep + LOO).
//
// CAMBIOS vs rsv_renewal.stan:
//   1. rho movido de `data` a `parameters` con constraint [0, 1)
//   2. Prior: Beta(alpha_rho, beta_rho), default Beta(2, 2).
//      Density -> 0 en {0, 1}, simetrico, mode 0.5, neutral.
//
// Compatible con MARSS (use_filter=1, sigma_obs como matriz P x T data)
// y con FSV-style (use_filter=0, sigma_obs_pl[p] como parametro).
//
// Ver rsv_renewal_freesigma_vec_rhoest.stan para la version FSV pura.
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
  // CAMBIO 1: rho ya NO esta en data
  vector<lower=1>[P] It_init;
  matrix<lower=0>[P, T] y;
  array[P, T] int<lower=0, upper=1> y_obs;
  real<lower=0> mean_gen_wk;
  real<lower=0> sd_gen_wk;
  int<lower=0, upper=1> use_filter;
  matrix<lower=0>[P, T] sigma_obs;
  // Hyper-parametros del prior Beta(alpha_rho, beta_rho).
  // Default Beta(2,2): alpha_rho = 2, beta_rho = 2 (simetrico, mode 0.5).
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
}

parameters {
  vector<lower=-20, upper=20>[P] log_beta;
  matrix<lower=-20, upper=20>[P, T] lam;
  matrix<lower=1e-6, upper=1e10>[P, T] It;
  vector<lower=0, upper=10>[P] sigma_obs_pl;
  // CAMBIO 1 (cont.): rho como parametro escalar en [0, 1)
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
  log_beta ~ normal(0, 4);
  sigma_obs_pl ~ normal(0, 2);  // half-Normal(0, 2^2) via <lower=0>; matches renewal_A_fsv_est.stan. Only active when use_filter==0.
  // CAMBIO 2: prior Beta(alpha_rho, beta_rho) para rho
  rho ~ beta(alpha_rho, beta_rho);

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

  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real s_eff = use_filter == 1 ? sigma_obs[p, t] : sigma_obs_pl[p];
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(s_eff);
        y[p, t] ~ lognormal(mu_obs, s_eff);
      }
    }
  }
}

generated quantities {
  matrix[P, T] log_lik;
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real s_eff = use_filter == 1 ? sigma_obs[p, t] : sigma_obs_pl[p];
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(s_eff);
        log_lik[p, t] = lognormal_lpdf(y[p, t] | mu_obs, s_eff);
      } else {
        log_lik[p, t] = 0;
      }
    }
  }
}
