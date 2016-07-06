functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  int xsz;
  vector[xsz] x;
}
model {}
generated quantities {
  int sz;
  int idx[6, 3];
  real ll;
  vector[p] v;
  matrix[p, p] Finv;
  matrix[m, p] K;
  vector[m] a;
  matrix[m, m] P;
  sz = ssm_filter_size(m, p);
  idx = ssm_filter_idx(m, p);
  ll = ssm_filter_get_loglik(x, m, p);
  v = ssm_filter_get_v(x, m, p);
  Finv = ssm_filter_get_Finv(x, m, p);
  K = ssm_filter_get_K(x, m, p);
  a = ssm_filter_get_a(x, m, p);
  P = ssm_filter_get_P(x, m, p);
}
