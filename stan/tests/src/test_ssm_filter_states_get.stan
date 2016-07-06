functions {
  #include ssm.stan
}
data {
  int m;
  int xsz;
  vector[xsz] x;
}
model {}
generated quantities {
  int sz;
  vector[m] a;
  matrix[m, m] P;
  sz = ssm_filter_states_size(m);
  a = ssm_filter_states_get_a(x, m);
  P = ssm_filter_states_get_P(x, m);
}
