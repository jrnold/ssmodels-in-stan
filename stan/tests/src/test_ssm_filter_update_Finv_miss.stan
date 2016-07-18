functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  int p_t;
  matrix[m, m] P;
  matrix[p, m] Z;
  matrix[p, p] H;
  int y_idx[p];
}
model{}
generated quantities {
  matrix[p, p] output;
  output = ssm_filter_update_Finv_miss(P, Z, H, p_t, y_idx);
}
