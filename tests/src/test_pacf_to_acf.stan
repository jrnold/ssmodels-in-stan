functions {
  #include ssm.stan
}
data {
  int n;
  vector[n] x;
}
model {}
generated quantities {
  vector[n] output;
  output = pacf_to_acf(x);
}
