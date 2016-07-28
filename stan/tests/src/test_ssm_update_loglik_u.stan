functions {
  #include ssm.stan
}
data {
  real v;
  real Finv;
}
model {}
generated quantities {
  real output;
  output = ssm_update_loglik_u(v, Finv);
}
