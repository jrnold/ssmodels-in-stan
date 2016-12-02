functions {
  #include ssm.stan
}
data {
  int n;
  matrix[n, n] A;
}
model {}
generated quantities {
  matrix[n, n] output;
  output = cholesky_decompose2(A);
}
