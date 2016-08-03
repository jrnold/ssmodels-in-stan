functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  real s;

  vector<lower = 0.0>[1] a1;
  cov_matrix[1] P1;
  real<lower = 0.0> y_scale;
}
transformed data {
  // system matrices
  matrix[1, 1] T[1];
  matrix[1, 1] Z[1];
  matrix[1, 1] R[1];
  vector[1] c[1];
  vector[1] d[1];
  int m;
  int p;
  int q;
  int filter_sz;
  m = 1;
  p = 1;
  q = 1;
  T[1, 1, 1] = 1.0;
  Z[1, 1, 1] = 1.0;
  R[1, 1, 1] = 1.0;
  c[1, 1] = 0.0;
  d[1, 1] = 0.0;
  filter_sz = ssm_filter_size(m, p);
}
parameters {
  real<lower = 0.0> sigma_eta;
  real<lower = 0.0> sigma_epsilon;
  vector<lower = 0.0>[n] lambda;
}
transformed parameters {
  matrix[1, 1] H[1];
  matrix[1, 1] Q[n];
  H[1, 1, 1] = pow(sigma_epsilon, 2);
  for (t in 1:n) {
    Q[t, 1, 1] = pow(sigma_epsilon * sigma_eta * lambda[t], 2);
  }
}
model {
  target += ssm_lp(y, d, Z, H, c, T, R, Q, a1, P1);
  sigma_epsilon ~ cauchy(0.0, y_scale);
  sigma_eta ~ cauchy(0.0, 1.0);
  lambda ~ cauchy(0.0, s);
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[1] alpha[n];
  // filtering
  filtered = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
  // sampling states
  alpha = ssm_simsmo_states_rng(filtered, d, Z, H, c, T, R, Q, a1, P1);
}
