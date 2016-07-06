functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  vector[p] v;
  matrix[p, p] Finv;
}
model {}
generated quantities {
  real output;
  output = ssm_filter_update_ll(v, Finv);
}
