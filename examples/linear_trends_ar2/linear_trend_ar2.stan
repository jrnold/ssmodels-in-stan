functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  vector[4] a1;
  cov_matrix[4] P1;
}
transformed data {
  matrix[1, 4] Z;
  matrix[4, 3] R;
  matrix[1, 1] H;
  vector[4] c;
  vector[1] d;
  d = rep_vector(0.0, 1);
  c = rep_vector(0.0, 4);
  R = rep_matrix(0.0, 4, 3);
  for (i in 1:3) {
    R[i, i] = 1.0;
  }
  Z[1, 1] = 1.0;
  Z[1, 2] = 0.0;
  Z[1, 3] = 1.0;
  Z[1, 4] = 0.0;
  H = rep_matrix(0.0, 1, 1);
}
parameters {
  vector<lower = 0.0>[3] sigma_eta;
  vector[2] phi_raw;
}
transformed parameters {
  vector[2] phi;
  phi = constrain_stationary(phi_raw);
}
model {
  matrix[4, 4] T;
  matrix[3, 3] Q;
  Q = rep_matrix(0.0, 3, 3);
  T = rep_matrix(0.0, 4, 4);
  T[1, 1] = 1.0;
  T[1, 2] = 1.0;
  T[2, 2] = 1.0;
  T[3, 3] = phi[1];
  T[4, 3] = phi[2];
  for (i in 1:3) {
    Q[i, i] = pow(sigma_eta[i], 2);
  }
  target += ssm_constant_lpdf(y | c, Z, H, d, T, R, Q, a1, P1);
}
