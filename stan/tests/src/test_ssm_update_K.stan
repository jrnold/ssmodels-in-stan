functions {
  #include ssm.stan
}
data {
  int m;
  int p;
  matrix[m, m] P;
  matrix[p, m] Z;
  matrix[m, m] T;
  matrix[p, p] Finv;
}
model {}
generated quantities {
  matrix[m, p] output;
  output = ssm_update_K(P, Z, T, Finv);
}
