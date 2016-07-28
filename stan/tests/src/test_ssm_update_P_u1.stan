functions {
  #include ssm.stan
}
data {
  int m;
  matrix[m, m] P;
  real Finv;
  vector[m] K;
}
model{}
generated quantities {
  matrix[m, m] output;
  output = ssm_update_P_u1(P, Finv, K);
}
