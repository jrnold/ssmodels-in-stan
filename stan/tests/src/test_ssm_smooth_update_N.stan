functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  matrix[m, m] N;
  matrix[p, m] Z;
  matrix[p, p] Finv;
  matrix[m, m] L;
}
model {}
generated quantities {
  matrix[m, m] output;
  output = ssm_smooth_update_N(N, Z, Finv, L);
}
