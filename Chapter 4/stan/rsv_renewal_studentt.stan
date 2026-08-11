// rsv_renewal_studentt.stan
// Variante MARSS-filtered (Specification B) con observation noise Student-t
// en escala log. Equivalente al rsv_renewal.stan original pero con la
// likelihood cambiada de lognormal a log-Student-t.
//
// MOTIVACIÓN:
//   Klaassen et al. (2024, Environ Res 240) modela wastewater data con
//   Student-t(df=10) para acomodar outliers de batch effects. Aquí adoptamos
//   el mismo principio en escala log, con sigma_obs[p,t] de MARSS Stage 1
//   como data (matriz P × T, no parámetro).
//
// CAMBIOS vs rsv_renewal.stan:
//   - nu_obs: nuevo data input (e.g., 10) — Student-t degrees of freedom
//   - likelihood:
//       Antes:  y ~ lognormal(log(psi) - 0.5 σ², σ)        donde σ = sigma_obs[p,t]
//       Ahora:  log(y) ~ student_t(nu, log(psi) - 0.5 σ², σ) - log(y)  (Jacobiano)
//
// REFERENCIA:
//   Klaassen, F. et al. (2024). "Predictive power of wastewater for
//   nowcasting infectious disease transmission." Environ Res 240(Pt 2): 117395.

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
  int<lower=0, upper=1> use_filter;          // 1 = use sigma_obs (MARSS), 0 = sigma_obs_pl (FSV-style)
  matrix<lower=0>[P, T] sigma_obs;           // P × T noise matrix from MARSS Stage 1
  real<lower=2> nu_obs;                       // NUEVO: Student-t degrees of freedom (e.g., 10)
}

parameters {
  vector<lower=-20, upper=20>[P] log_beta;
  matrix<lower=-20, upper=20>[P, T] lam;
  matrix<lower=1e-6, upper=1e10>[P, T] It;
  vector<lower=0, upper=10>[P] sigma_obs_pl;  // used only when use_filter=0
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
  // ---- PRIORS ---------------------------------------------------------------
  log_beta ~ normal(0, 4);
  sigma_obs_pl ~ lognormal(log(0.5), sqrt(10));   // only used if use_filter == 0

  for (p in 1:P) lam[p, 1] ~ normal(0, 0.5);
  for (t in 2:T) {
    for (p in 1:P) {
      real mu_lam = 0;
      for (j in 1:P) mu_lam += W[p, j] * lam[j, t - 1];
      mu_lam = (1 - rho) * lam[p, t - 1] + rho * mu_lam;
      lam[p, t] ~ normal(mu_lam, 0.5);
    }
  }

  // ---- LATENT It state (Poisson → LogNormal moment-matched relaxation) -----
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

  // ---- OBSERVATION MODEL — CAMBIO PRINCIPAL: Student-t en log(y) ----------
  // Y sigue log-Student-t:
  //   log p_Y(y) = student_t_lpdf(log(y) | nu, mu, sigma) - log(y)
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real s_eff  = use_filter == 1 ? sigma_obs[p, t] : sigma_obs_pl[p];
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(s_eff);
        target += student_t_lpdf(log(y[p, t]) | nu_obs, mu_obs, s_eff)
                - log(y[p, t]);   // Jacobiano del cambio z = log(y)
      }
    }
  }
}

generated quantities {
  matrix[P, T] log_lik;
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real s_eff  = use_filter == 1 ? sigma_obs[p, t] : sigma_obs_pl[p];
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(s_eff);
        // log_lik en escala Y (incluye Jacobiano), comparable directamente
        // con el modelo lognormal vía PSIS-LOO.
        log_lik[p, t] = student_t_lpdf(log(y[p, t]) | nu_obs, mu_obs, s_eff)
                      - log(y[p, t]);
      } else {
        log_lik[p, t] = 0;
      }
    }
  }
}
