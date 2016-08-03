functions {
  #include ssm.stan
}
data {
  int n;
  int p;
  int m;
  int q;
  vector[p] y[n];
  int d_sz;
  vector[p] d[d_sz];
  int Z_sz;
  matrix[p, m] Z[Z_sz];
  int H_sz;
  matrix[p, p] H[H_sz];
  int c_sz;
  vector[p] c[c_sz];
  int T_sz;
  matrix[m, m] T[T_sz];
  int R_sz;
  matrix[m, q] R[R_sz];
  int Q_sz;
  matrix[q, q] Q[Q_sz];
  vector[m] a1;
  matrix[m, m] P1;
}
model {
}
generated quantities {
  real output;
  output = ssm_loglik(y, d, Z, H, c, T, R, Q, a1, P1);
}
