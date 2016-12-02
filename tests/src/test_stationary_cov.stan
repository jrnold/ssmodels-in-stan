functions {
  #include ssm.stan
}
data {
  int m;
  matrix[m, m] T;
  matrix[m, m] RQR;
}
model {}
generated quantities {
  matrix[m, m] output;
  output = stationary_cov(T, RQR);
}
