functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  int q;
  int<lower = 1> n;
  vector[p] y[n];
  vector[p] d[1];
  matrix[p, m] Z[1];
  matrix[p, p] H[1];
  vector[m] c[1];
  matrix[m, m] T[1];
  matrix[m, q] R[1];
  matrix[q, q] Q[1];
  vector<lower = 0.0>[m] a1;
  cov_matrix[m] P1;
}
transformed data {
  int output_size;
  output_size = ssm_filter_size(m, p);
}
model {}
generated quantities {
  vector[output_size] output[n];
  output = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
}
