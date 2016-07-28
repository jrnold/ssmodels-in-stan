functions {
  #include ssm.stan
}
data {
  int m;
  vector[m] a;
  vector[m] c;
  matrix[m, m] T;
}
model{}
generated quantities {
  vector[m] output;
  output = ssm_update_a_u2(a, c, T);
}
