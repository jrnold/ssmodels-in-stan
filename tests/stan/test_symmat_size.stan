functions {
  #include ssm.stan
}
data {
  int n;
}
model {}
generated quantities {
  int output;
  output = symmat_size(n);
}
