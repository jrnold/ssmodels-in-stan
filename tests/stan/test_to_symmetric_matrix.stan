functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  matrix[m, n] input;
}
model {}
generated quantities {
  matrix[m, m] output;
  output = to_symmetric_matrix(input);
}
