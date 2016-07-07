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
  output = acf_to_pacf(x);
}
