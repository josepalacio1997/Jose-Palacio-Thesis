// rsv_renewal_freesigma_vec_studentt.stan
// Variante de rsv_renewal_freesigma_vec.stan con observation noise
// Student-t en escala log (i.e., Y sigue una log-Student-t).
//
// MOTIVACIÓN:
//   Klaassen et al. (2024, Environ Res 240) modela wastewater data con
//   Student-t(df=10) en escala lineal. Aquí adoptamos el mismo principio
//   en escala log para mantener consistencia con la parametrización
//   lognormal del modelo original. Las colas pesadas (vs Lognormal)
//   acomodan outliers de batch effects sin re-fittear con datos limpios.
//
// CAMBIOS vs rsv_renewal_freesigma_vec.stan:
//   - sigma_obs_pl: vector[P] sin cambio (per-planta, libre)
//   - nu_obs: nuevo parámetro fijo (data input) = 10 por default
//   - likelihood:
//       Antes:  y ~ lognormal(log(psi) - 0.5 sigma^2, sigma)
//       Ahora:  log(y) ~ student_t(nu_obs, log(psi) - 0.5 sigma^2, sigma)
//                + Jacobian (-log y) para la comparación LOO con Lognormal
//
// REFERENCIA:
//   Klaassen, F. et al. (2024). "Predictive power of wastewater for
//   nowcasting infectious disease transmission: a retrospective case
//   study of five sewershed areas in Louisville, Kentucky."
//   Environmental Research 240(Pt 2): 117395.

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
  real<lower=2> nu_obs;            // NUEVO: grados de libertad Student-t (e.g., 10)
}

parameters {
  vector<lower=-20, upper=20>[P] log_beta;
  matrix<lower=-20, upper=20>[P, T] lam;
  matrix<lower=1e-6, upper=1e10>[P, T] It;
  vector<lower=0, upper=10>[P] sigma_obs_pl;   // sin cambio vs FSV original
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
    for (t in 2:T) {
      Rt[p, t] = euler_lotka_rt(lam[p, t], mean_gen_wk, sd_gen_wk);
      real conv = 0;
      for (s in 1:min(G, t - 1)) {
        conv += g[s] * It[p, t - s];
      }
      iota[p, t] = Rt[p, t] * fmax(conv, 1e-10);
    }
  }

  // Mixing espacial
  for (t in 1:T) {
    for (p in 1:P) {
      real local_sum = 0;
      for (q in 1:P) local_sum += W[p, q] * It[q, t];
      I_mix[p, t] = (1 - rho) * It[p, t] + rho * local_sum;
    }
  }

  // Convolución con shedding kernel
  for (p in 1:P) {
    for (t in 1:T) {
      real shed = 0;
      for (k in 0:min(D - 1, t - 1)) {
        shed += d[k + 1] * I_mix[p, t - k];
      }
      psi[p, t] = beta[p] * shed;
    }
  }
}

model {
  // Priors
  log_beta ~ normal(0, 4);
  sigma_obs_pl ~ normal(0, 2);   // Half-Normal por componentwise

  // Latent AR-like + prior sobre It[1]
  for (p in 1:P) {
    real lam_p = fmin(fmax(It_init[p], 10.0), 1e8);
    target += lognormal_lpdf(It[p, 1] | log(lam_p) - 0.5/lam_p,
                                         1/sqrt(lam_p));
  }
  for (p in 1:P) for (t in 1:T) {
    lam[p, t] ~ normal(0.1, 0.5);
  }
  for (p in 1:P) for (t in 2:T) {
    real lam_p = fmin(fmax(iota[p, t], 10.0), 1e8);
    target += lognormal_lpdf(It[p, t] | log(lam_p) - 0.5/lam_p,
                                         1/sqrt(lam_p));
  }

  // OBSERVATION MODEL — CAMBIO PRINCIPAL: Student-t en log(y) con Jacobian
  // y ~ log-Student-t(nu_obs, log(psi) - 0.5*sigma^2, sigma_obs_pl[p])
  for (p in 1:P) {
    for (t in 1:T) {
      if (y_obs[p, t] == 1) {
        real mu_obs = log(psi[p, t] + 1e-6) - 0.5 * square(sigma_obs_pl[p]);
        // Equivalente a "Y sigue log-Student-t":
        //   log p_Y(y) = log p_Z(log(y)) - log(y)
        //             = student_t_lpdf(log(y) | nu, mu, sigma) - log(y)
        target += student_t_lpdf(log(y[p, t]) | nu_obs, mu_obs, sigma_obs_pl[p])
                - log(y[p, t]);
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
        // log_lik en escala Y (comparable directamente con lognormal model
        // para LOO):
        log_lik[p, t] = student_t_lpdf(log(y[p, t]) | nu_obs, mu_obs,
                                       sigma_obs_pl[p])
                      - log(y[p, t]);
      } else {
        log_lik[p, t] = 0;
      }
    }
  }
}
