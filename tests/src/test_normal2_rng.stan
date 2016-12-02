functions {
  #include ssm.stan
}
data {
  int n;
  real mu;
  real sigma;
}
model {}
generated quantities {
  vector[n] output;
  for (i in 1:n) {
    output[i] = normal2_rng(mu, sigma);
  }
}
