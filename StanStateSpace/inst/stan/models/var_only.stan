functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  int<lower = 1> m;
  int<lower = 1> q;
  vector[1] y[n];
  vector[m] a1;
  cov_matrix[m] P1;
  matrix[m, m] T;
  matrix[1, m] Z;
  matrix[m, q] R;
  vector[m] c;
  vector[1] d;
}
transformed data {
  int filter_sz;
  int p;
  real y_sd;
  p = 1;
  filter_sz = ssm_filter_size(m, p);
  {
    vector[n] yvec;
    for (i in 1:n) {
      yvec[i] = y[i][1];
    }
    y_sd = sd(yvec);
  }
}
parameters {
  vector<lower = 0.0>[q] sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
  matrix[1, 1] H;
  matrix[q, q] Q;
  H = rep_matrix(pow(sigma_epsilon, 2), 1, 1);
  Q = rep_matrix(0.0, q, q);
  for (i in 1:q) {
    Q[i, i] = pow(sigma_eta[i] * sigma_epsilon, 2);
  }
}
model {
  target += ssm_constant_lp(d, Z, H, c, T, R, Q, a1, P1);
  sigma_epsilon ~ cauchy(0.0, y_sd);
  sigma_eta ~ cauchy(0.0, 1.0);
}
generated quantities {
  vector[m] alpha[n];
  vector[m] alpha_hat[n];
  vector[filter_sz] filtered[n];
  filtered = ssm_filter_miss(y, d, Z, H, c, T, R, Q, a1, P1, p_t, y_idx);
  alpha_hat = ssm_smooth_states_mean(filtered, Z, c, T, R, Q);
  alpha = ssm_simsmo_states_rng(filtered, d, Z, H, c, T, R, Q, a1, P1);
}
