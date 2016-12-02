functions {
  #include ssm.stan
}
data {
  int n;
}
model {}
generated quantities {
  int output;
  output = find_symmat_dim(n);
}
