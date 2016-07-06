functions {
  #include ssm.stan
}
data {
  int p;
  int xsz;
  vector[xsz] x;
}
model {}
generated quantities {
  int sz;
  vector[p] eps_mean;
  matrix[p, p] eps_var;
  sz = ssm_smooth_eps_size(p);
  eps_mean = ssm_smooth_eps_get_mean(x, p);
  eps_var = ssm_smooth_eps_get_var(x, p);
}
