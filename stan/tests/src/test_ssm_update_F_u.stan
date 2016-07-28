functions {
  #include ssm.stan
}
data {
  int m;
  matrix[m, m] P;
  row_vector[m] Z;
  real H;
}
model{}
generated quantities {
  real output;
  output = ssm_update_F_u(P, Z, H);
}
