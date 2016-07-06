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
  vector[m] alpha;
  matrix[m, m] V;
  sz = ssm_smooth_state_size(m);
  alpha = ssm_smooth_state_get_mean(x, m);
  V = ssm_smooth_state_get_var(x, m);
}
