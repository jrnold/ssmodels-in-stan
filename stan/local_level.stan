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
  T[1, 1] <- 1.0;
  Z[1, 1] <- 1.0;
  R[1, 1] <- 1.0;
  c[1] <- 0.0;
  d[1] <- 0.0;
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
    H <- rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    Q <- rep_matrix(pow(sigma_eta, 2), 1, 1);
    ssm_lp(y, c, Z, H, d, T, R, Q, a1, P1);
  }
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[m + m * m] eta[n];
  vector[1] eps[n];
  vector[m + m * m] alpha[n];
  vector[1] alpha2[n];
  {
    matrix[1, 1] H;
    matrix[1, 1] Q;
    H <- rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    Q <- rep_matrix(pow(sigma_eta, 2), 1, 1);
    filtered <- ssm_filter(y, c, Z, H, d, T, R, Q, a1, P1);
    eps <- ssm_smooth_eps(filtered, H, Z, T);
    eta <- ssm_smooth_eta(filtered, Z, T, R, Q);
    alpha <- ssm_smooth_state(filtered, Z, T);
    alpha2 <- ssm_smooth_faststate(filtered, Z, T, R, Q);
  }
}
