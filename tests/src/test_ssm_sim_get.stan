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
  int idx[2, 3];
  vector[p] y;
  vector[m] a;
  sz = ssm_sim_size(m, p, q);
  idx = ssm_sim_idx(m, p, q);
  y = ssm_sim_get_y(x, m, p, q);
  a = ssm_sim_get_a(x, m, p, q);
}
