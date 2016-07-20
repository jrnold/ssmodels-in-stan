functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  matrix[m, m] P;
  matrix[p, m] Z;
  matrix[p, p] H;
}
model {}
generated quantities {
  matrix[p, p] output;
  output = ssm_update_F(P, Z, H);
}
