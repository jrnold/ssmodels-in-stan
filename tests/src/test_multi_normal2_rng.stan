functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  vector[m] mu;
  matrix[m, m] Sigma;
}
model {}
generated quantities {
  vector[m] output[n];
  for (i in 1:n) {
    output[i] = multi_normal2_rng(mu, Sigma);
  }
}
