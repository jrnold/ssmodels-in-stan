functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  vector[m] a;
  vector[m] c;
  matrix[m, m] T;
  vector[p] v;
  matrix[m, p] K;
}
model {}
generated quantities {
  vector[m] output;
  output = ssm_update_a(a, c, T, v, K);
}
