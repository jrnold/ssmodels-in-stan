functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  vector[m] mu;
  matrix[m, m] L;
}
model {}
generated quantities {
  vector[m] output[n];
  for (i in 1:n) {
    output[i] = multi_normal_cholesky2_rng(mu, L);
  }
}
