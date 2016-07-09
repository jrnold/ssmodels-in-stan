functions {
  #include ssm.stan
}
data {
  int n;
  int m;
  int x[m];
}
model {}
generated quantities {
  int output[n];
  output = mask_indexes(x, n);
}
