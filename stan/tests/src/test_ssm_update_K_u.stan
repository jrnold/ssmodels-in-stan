functions {
  #include ssm.stan
}
data {
  int m;
  matrix[m, m] P;
  row_vector[m] Z;
  real Finv;
}
model{}
generated quantities {
  vector[m] Finv;
  output = ssm_update_K_u(P, Z, Finv);
}
