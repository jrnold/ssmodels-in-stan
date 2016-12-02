functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  vector[n * m] input;
}
model {}
generated quantities {
  matrix[m, n] output;
  output = to_matrix_colwise(input, m, n);
}
