functions {
  #include ssm.stan
}
data {
  int m;
  int n;
  int p;
  int q;
  matrix[p, q] x;
  int i[p];
  int j[q];
  real a;
}
model {}
generated quantities {
  matrix[m, n] output;
  output = fill_matrix(x, m, n, i, j, a);
}
