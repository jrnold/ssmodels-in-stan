functions {
  #include ssm.stan
}
data {
  int m;
  vector[m] a;
  row_vector[m] Z;
  real d;
  real y;
}
model{}
generated quantities {
  real output;
  output = ssm_update_v_u(y, a, d, Z);
}
