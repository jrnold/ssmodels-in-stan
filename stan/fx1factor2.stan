data {
  // dimensions
  int n; // number of observations
  int r; // number of variables
  int p; // number of states
  // observations
  matrix[r, n] y;
  // system matrices
  // observation equation
  matrix[r, p] F;
  // system equation
  matrix[p, p] G;
  // initial conditions
  vector[p] m0;
  cov_matrix[p] C0;
}
parameters {
  cov_matrix[p] W;
  vector<lower = 0.0>[r] V;
}
model {
  y ~ gaussian_dlm_obs(F ', G, V, W, m0, C0);
}
