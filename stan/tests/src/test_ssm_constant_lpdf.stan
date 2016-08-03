functions {
  #include ssm.stan
}
data {
  int n;
  int p;
  int m;
  int q;
  vector[p] y[n];
  vector[p] d;
  matrix[p, m] Z;
  matrix[p, p] H;
  vector[p] c;
  matrix[m, m] T;
  matrix[m, q] R;
  matrix[q, q] Q;
  vector[m] a1;
  matrix[m, m] P1;
}
parameters {
}
transformed parameters {
}
model {
}
generated quantities {
  real output;
  output = ssm_constant_loglik(y, d, Z, H, c, T, R, Q, a1, P1);
}
