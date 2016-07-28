functions {
  #include ssm.stan
}
data {
  int m;
  matrix[m, m] P;
  matrix[m, m] T;
  matrix[m, m] RQR;
}
model{}
generated quantities {
  matrix[m, m] output;
  output = ssm_update_P_u2(P, T, RQR);
}
