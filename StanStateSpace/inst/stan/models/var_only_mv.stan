functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  int<lower = 1> m;
  int<lower = 1> p;
  int<lower = 1, upper = m> q;
  vector[p] y[n];

  // State Space
  int<lower = 1, upper = n> d_sz;
  vector[p] d[d_sz];
  int<lower = 1, upper = n> Z_sz;
  matrix[p, m] Z[Z_sz];
  int<lower = 1, upper = n> c_sz;
  vector[m] c[d_sz];
  int<lower = 1, upper = n> T_sz;
  matrix[m, m] T[T_sz];
  int<lower = 1, upper = n> R_sz;
  matrix[m, q] R[R_sz];
  vector<lower = 0.>[m] a1;
  matrix[m, m] P1;
  // Hyperparameters
  real<lower=1.> Rho_eta_prior;
  real<lower=1.> Rho_epsilon_prior;
  vector<lower=0.>[q] tau_eta_prior;
  vector<lower=0.>[p] tau_epsilon_prior;
}
transformed data {
  int filter_sz;
  filter_sz = ssm_filter_size(m, p);
}
parameters {
  vector<lower = 0.>[p] tau_epsilon;
  vector<lower = 0.>[q] tau_eta;
  corr_matrix[p] Rho_epsilon;
  corr_matrix[q] Rho_eta;
}
transformed parameters {
  matrix[p, p] H[1];
  matrix[q, q] Q[1];
  H[1] = quad_form_diag(Rho_epsilon, tau_epsilon);
  Q[1] = quad_form_diag(Rho_eta, tau_eta .* tau_epsilon);
}
model {
  target += ssm_lpdf(y | d, Z, H, c, T, R, Q, a1, P1);
  Rho_epsilon ~ lkj_corr(Rho_epsilon_prior);
  Rho_eta ~ lkj_corr(Rho_eta_prior);
  tau_epsilon ~ cauchy(0., tau_epsilon_prior);
  tau_eta ~ cauchy(0., tau_eta_prior);
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[m] alpha_mean[n];
  vector[m] alpha[n];
  {
    // filtering
    filtered = ssm_filter(y, d, Z, H,
                          c, T, R, Q,
                          a1, P1);
    // mean states
    alpha_mean = ssm_smooth_states_mean(filtered, Z, c, T, R, Q);
    // sampling states
    alpha = ssm_simsmo_states_rng(filtered, d, Z, H,
                                  c, T, R, Q, a1, P1);
  }
}
