functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  matrix[m, m] P;
  matrix[p, p] Finv;
  matrix[p, m] Z;
}
model {}
generated quantities {
  matrix[m, m] output;
  output = ssm_filter_states_update_P(P, Z, Finv);
}
