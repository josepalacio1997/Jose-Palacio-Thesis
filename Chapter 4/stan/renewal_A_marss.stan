
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
  array[P, T] int<lower=0, upper=1> y_obs;   // new Stan array syntax (cmdstan >= 2.32)
  real<lower=0> mean_gen_wk;
  real<lower=0> sd_gen_wk;
  int<lower=0, upper=1> use_filter;
  matrix<lower=0>[P, T] sigma_obs;
}

parameters {
  vector<lower=-20, upper=20>[P] log_beta;          // beta in [exp(-10), exp(10)]
  matrix<lower=-20, upper=20>[P, T] lam;            // keeps Rt bounded
  matrix<lower=1e-6, upper=1e10>[P, T] It;          // cap from above to avoid It -> inf
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
  sigma_obs_pl ~ normal(0, 2);  // half-Normal(0, 2^2) via <lower=0>; matches renewal_A_fsv.stan. Only active when use_filter==0.

  for (p in 1:P) lam[p, 1] ~ normal(0, 0.5);
  for (t in 2:T) {
    for (p in 1:P) {
      real mu_lam = 0;
      for (j in 1:P) mu_lam += W[p, j] * lam[j, t - 1];
      mu_lam = (1 - rho) * lam[p, t - 1] + rho * mu_lam;
      lam[p, t] ~ normal(mu_lam, 0.5);
    }
  }

  // Poisson -> LogNormal moment-matched relaxation. We floor lam_p at 1.0
  // because mu_ln = log(lam) - 0.5/lam blows up to -inf as lam -> 0, which
  // NUTS hits during warmup when It (lower bound 1e-6) gets proposed near
  // zero. A floor of 1.0 keeps the relaxation numerically well-defined
  // without distorting the prior at biologically realistic incidence levels.
  // Poisson -> LogNormal moment-matched relaxation, clamped at BOTH ends.
  // Floor at 1.0 prevents log(lam) - 0.5/lam from going to -inf when iota -> 0.
  // Cap at 1e8 prevents log(lam) from going to +inf when iota -> +inf during
  // wild NUTS warmup proposals. Both blows are observed in practice and were
  // producing "Location parameter is inf" rejections at line 121.
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

