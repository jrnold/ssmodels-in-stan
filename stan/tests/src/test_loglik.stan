functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  int q;
  int<lower = 1> n;
  vector[p] y[n];
  vector[p] d;
  matrix[p, m] Z;
  matrix[p, p] H;
  vector[m] c;
  matrix[m, m] T;
  matrix[m, q] R;
  matrix[q, q] Q;
  vector<lower = 0.0>[m] a1;
  cov_matrix[m] P1;
}
generated quantities {
  real ll;
  ll = ssm_constant_lpdf(y, d, Z, H, c, T, R, Q, a1, P1);
}
