functions {
  #include ssm.stan
}
data {
  int m;
  row_vector[m] Z;
  vector[m] K;
}
model{}
generated quantities {
  matrix[m, m] output;
  output = ssm_update_L_u(Z, K);
}
