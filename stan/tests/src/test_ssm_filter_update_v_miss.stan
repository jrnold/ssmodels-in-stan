functions {
  #include ssm.stan
}
data {
  int p;
  int p_t;
  int m;
  vector[p] y;
  vector[m] a;
  vector[p] d;
  matrix[p, m] Z;
  int y_idx[p];
}
model {}
generated quantities {
  vector[p] output;
  output = ssm_filter_update_v_miss(y, a, d, Z, p_t, y_idx);
}
