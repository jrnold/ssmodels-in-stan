functions {
  #include ssm.stan
}
data {
  int p;
  int p_t;
  int m;
  vector[p] v;
  matrix[p, p] Finv;
  int y_idx[p];
}
model {}
generated quantities {
  real output;
  output = ssm_update_ll_miss(v, Finv, p_t, y_idx);
}
