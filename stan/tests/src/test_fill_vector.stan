functions {
  #include ssm.stan
}
data {
  int m;
  int n;
  vector[m] x;
  int i[m];
  real a;
}
model {}
generated quantities {
  vector[n] output;
  output = fill_vector(x, n, i, a);
}
