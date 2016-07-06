functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  matrix[m, m] P;
  matrix[p, m] Z;
  matrix[m, m] T;
  matrix[m, m] RQR;
  matrix[m, p] K;
}
model {
}
generated quantities {
  matrix[m, m] output;
  output = ssm_filter_update_P(P, Z, T, RQR, K);
}
