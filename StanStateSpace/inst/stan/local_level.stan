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
  real y_sd;
  m = 1;
  p = 1;
  q = 1;
  T[1, 1] = 1.0;
  Z[1, 1] = 1.0;
  R[1, 1] = 1.0;
  c[1] = 0.0;
  d[1] = 0.0;
  {
    vector[n] yvec;
    for (i in 1:n) {
      yvec[i] = y[i][1];
    }
    y_sd = sd(yvec);
  }
}
parameters {
  real<lower = 0.0> sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
  matrix[1, 1] H;
  matrix[1, 1] Q;
  H = rep_matrix(pow(sigma_epsilon, 2), 1, 1);
  Q = rep_matrix(pow(sigma_eta * sigma_epsilon, 2), 1, 1);
}
model {
  y ~ ssm_constant_lpdf(d, Z, H, c, T, R, Q, a1, P1);
  sigma_epsilon ~ cauchy(0.0, y_sd);
  sigma_eta ~ cauchy(0.0, 1.0);
}
generated quantities {
  vector[6] filtered[n];
  vector[2] eta_hat[n];
  vector[2] eps_hat[n];
  vector[2] alpha_hat[n];
  vector[1] alpha_hat_fast[n];
  {
    // Filtered data
    filtered = ssm_filter(y,
                          rep_array(d, 1), rep_array(Z, 1), rep_array(H, 1),
                          rep_array(c, 1), rep_array(T, 1), rep_array(R, 1),
                          rep_array(Q, 1), a1, P1);
    // Smoothed values
    // state disturbances (full)
    alpha_hat = ssm_smooth_state(filtered, rep_array(Z, 1), rep_array(T, 1));
    // state distrurbances calculated using fast smoother
    alpha_hat_fast = ssm_smooth_faststate(filtered, rep_array(c, 1), rep_array(Z, 1),
                                  rep_array(T, 1), rep_array(R, 1),
                                  rep_array(Q, 1));
    // observation disturbance
    eps_hat = ssm_smooth_eps(filtered,
                             rep_array(H, 1), rep_array(Z, 1),
                             rep_array(T, 1));
    // state disturbances
    eta_hat = ssm_smooth_eta(filtered, rep_array(Z, 1), rep_array(T, 1),
                             rep_array(R, 1), rep_array(Q, 1));
    // Simulations

  }
}
