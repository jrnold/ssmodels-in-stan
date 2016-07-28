functions {
  #include ssm.stan
}
data {
  int m;
  vector[m] a;
  vector[m] K;
  real v;
}
model{}
generated quantities {
  vector[m] output;
  output = ssm_update_a_u1(a, v, K);
}
