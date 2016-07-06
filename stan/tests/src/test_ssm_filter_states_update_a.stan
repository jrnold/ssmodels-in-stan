functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  vector[m] a;
  matrix[m, m] P;
  matrix[p, p] Finv;
  vector[p] v;
  matrix[p, m] Z;
}
model {}
generated quantities {
  vector[m] output;
  output = ssm_filter_states_update_a(a, P, Z, v, Finv);
}
