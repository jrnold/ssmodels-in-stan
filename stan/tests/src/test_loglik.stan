functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  matrix[1, 1] H;
  matrix[1, 1] Q;
  vector<lower = 0.0>[1] a1;
  cov_matrix[1] P1;
}
transformed data {
  // system matrices
  matrix[1, 1] T;
  matrix[1, 1] Z;
  matrix[1, 1] R;
  vector[1] c;
  vector[1] d;
  matrix[1, 1] H;
  matrix[1, 1] Q;
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
}
transformed parameters {
}
model {
}
generated quantities {
  real ll;
  ll = ssm_constant_lpdf(y, d, Z, H, c, T, R, Q, a1, P1);
}
