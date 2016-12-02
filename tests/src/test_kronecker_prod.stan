functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  int p;
  int q;
  matrix[m, n] A;
  matrix[p, q] B;
}
model {}
generated quantities {
  matrix[m * p, n * q] output;
  output = kronecker_prod(A, B);
}
