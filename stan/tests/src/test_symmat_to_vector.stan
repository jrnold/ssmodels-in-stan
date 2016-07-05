functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  matrix[m, n] x;
}
transformed data {
  int nm;
  nm = symmat_size(min(n, m));
}
model {}
generated quantities {
  vector[nm] output;
  output = symmat_to_vector(x);
}
