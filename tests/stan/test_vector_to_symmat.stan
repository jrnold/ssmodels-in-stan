functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  vector[m] x;
}
model {}
generated quantities {
  matrix[n, n] output;
  output = vector_to_symmat(x, n);
}
