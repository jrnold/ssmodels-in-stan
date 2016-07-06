functions {
  #include ssm.stan
}
data {
  int m;
  int n;
  matrix[m, n] A;
  matrix[m, n] B;
}
model {}
generated quantities {
  real output;
  output = matrix_diff(A, B);
}
