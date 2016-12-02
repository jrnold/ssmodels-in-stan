functions {
  #include ssm.stan
}
data {
  real x;
  int m;
  int n;
  int k;
}
model {}
generated quantities {
  matrix[m, n] output;
  output = rep_diagonal_matrix(x, m, n, k);
}
