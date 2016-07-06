functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  vector[m] r;
  matrix[p, m] Z;
  vector[p] v;
  matrix[p, p] Finv;
  matrix[m, m] L;
}
model {}
generated quantities {
  vector[m] output;
  output = ssm_smooth_update_r(r, Z, v, Finv, L);
}
