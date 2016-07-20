functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  matrix[p, m] Z;
  matrix[m, m] T;
  matrix[m, p] K;
}
model {}
generated quantities {
  matrix[m, m] output;
  output = ssm_update_L(Z, T, K);
}
