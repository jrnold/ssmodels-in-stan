functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  vector<lower = 0.0>[1] a1;
  cov_matrix[1] P1;
}
transformed data {
  matrix[1, 1] T;
  matrix[1, 1] Z;
  matrix[1, 1] R;
  vector[1] c;
  vector[1] d;
  int m;
  int p;
  int q;
  m = 1;
  p = 1;
  q = 1;
  T[1, 1] = 1.0;
  Z[1, 1] = 1.0;
  R[1, 1] = 1.0;
  c[1] = 0.0;
  d[1] = 0.0;
}
parameters {
  real<lower = 0.0> sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
}
model {
  {
    matrix[1, 1] H;
    matrix[1, 1] Q;
    H = rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    Q = rep_matrix(pow(sigma_eta, 2), 1, 1);
    target += ssm_constant_lpdf(y | c, Z, H, d, T, R, Q, a1, P1);
  }
}
generated quantities {
  vector[6] filtered[n];
  vector[2] eta[n];
  vector[2] eps[n];
  vector[2] alpha[n];
  vector[1] alpha2[n];
  {
    matrix[1, 1] H;
    matrix[1, 1] Q;
    H = rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    Q = rep_matrix(pow(sigma_eta, 2), 1, 1);
    filtered = ssm_filter(y, rep_array(d, 1), rep_array(Z, 1), rep_array(H, 1),
                          rep_array(c, 1), rep_array(T, 1), rep_array(R, 1), rep_array(Q, 1), a1, P1);
    eps = ssm_smooth_eps(filtered, rep_array(H, 1), rep_array(Z, 1), rep_array(T, 1));
    eta = ssm_smooth_eta(filtered, rep_array(Z, 1), rep_array(T, 1), rep_array(R, 1), rep_array(Q, 1));
    alpha = ssm_smooth_state(filtered, rep_array(Z, 1), rep_array(T, 1));
    alpha2 = ssm_smooth_faststate(filtered, rep_array(c, 1), rep_array(Z, 1), rep_array(T, 1),
                                  rep_array(R, 1), rep_array(Q, 1));
  }
}
