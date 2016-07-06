functions {
  #include ssm.stan
}
data {
  int q;
  int xsz;
  vector[xsz] x;
}
model {}
generated quantities {
  int sz;
  vector[q] eta_mean;
  matrix[q, q] eta_var;
  sz = ssm_smooth_eta_size(q);
  eta_mean = ssm_smooth_eta_get_mean(x, q);
  eta_var = ssm_smooth_eta_get_var(x, q);
}
