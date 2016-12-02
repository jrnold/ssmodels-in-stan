functions {
  #include ssm.stan
}
data {
  int m;
  int n;
  matrix[m, m] A;
}
model {}
generated quantities {
  matrix[m, m] output;
  output = matrix_pow(A, n);
}
