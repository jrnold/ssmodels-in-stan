functions {
  #include ssm.stan
}
data {
  int p;
  int m;
  vector[p] y;
  vector[m] a;
  vector[p] d;
  matrix[p, m] Z;
}
model {}
generated quantities {
  vector[p] output;
  output = ssm_filter_update_v(y, a, d, Z);
}
