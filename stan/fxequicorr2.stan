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
transformed data {
}
parameters {
  vector<lower = 0.0>[p] wdiag;
  vector<lower = 0.0>[r] vdiag;
  real<lower = -1, upper = 1> rho;
}
transformed parameters {
  // system variance matrix
  matrix[p, p] W;
  matrix[r, r] V;
  W <- diag_matrix(wdiag);
  for (i in 1:r) {
    for (j in 1:r) {
      if (i == j) {
	V[i, i] <- pow(vdiag[i], 2);
      } else {
	V[i, j] <- vdiag[i] * vdiag[j] * rho;
      }
    }
  }
}
model {
  y ~ gaussian_dlm_obs(F', G, V, W, m0, C0);
}
generated quantities {
}