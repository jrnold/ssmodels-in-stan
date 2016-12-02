functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  vector<lower = 0.0>[1] a1;
  cov_matrix[1] P1;
  real<lower = 0.0> y_scale;
}
transformed data {
  // system matrices
  matrix[1, 1] T;
  matrix[1, 1] Z;
  matrix[1, 1] R;
  vector[1] c;
  vector[1] d;
  int m;
  int p;
  int q;
  int filter_sz;
  m = 1;
  p = 1;
  q = 1;
  T[1, 1] = 1.0;
  Z[1, 1] = 1.0;
  R[1, 1] = 1.0;
  c[1] = 0.0;
  d[1] = 0.0;
  filter_sz = ssm_filter_size(m, p);
}
parameters {
  real<lower = 0.0> sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
  matrix[1, 1] H;
  matrix[1, 1] Q;
  H[1, 1] = pow(sigma_epsilon, 2);
  Q[1, 1] = pow(sigma_eta * sigma_epsilon, 2);
}
model {
  target += ssm_constant_lpdf(y | d, Z, H, c, T, R, Q, a1, P1);
  sigma_epsilon ~ cauchy(0.0, y_scale);
  sigma_eta ~ cauchy(0.0, 1.0);
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[1] alpha[n];
  vector[1] alpha_mean[n];
  // filtering
  filtered = ssm_filter(y,
    rep_array(d, 1),
    rep_array(Z, 1),
    rep_array(H, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1), a1, P1);
  alpha_mean = ssm_smooth_states_mean(filtered,
    rep_array(Z, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1));
  // sampling states
  alpha = ssm_simsmo_states_rng(filtered,
    rep_array(d, 1),
    rep_array(Z, 1),
    rep_array(H, 1),
    rep_array(c, 1),
    rep_array(T, 1),
    rep_array(R, 1),
    rep_array(Q, 1), a1, P1);
}
