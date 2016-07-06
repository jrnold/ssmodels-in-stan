functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  int q;
  int xsz;
  vector[xsz] x;
}
model {}
generated quantities {
  int sz;
  int idx[4, 3];
  vector[p] y;
  vector[m] a;
  vector[p] eps;
  vector[q] eta;
  sz = ssm_sim_size(m, p, q);
  idx = ssm_sim_idx(m, p, q);
  y = ssm_sim_get_y(x, m, p, q);
  a = ssm_sim_get_a(x, m, p, q);
  eps = ssm_sim_get_eps(x, m, p, q);
  eta = ssm_sim_get_eta(x, m, p, q);
}
