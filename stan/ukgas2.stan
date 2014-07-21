data {
  // dimensions
  int n; // number of observations
  int r; // number of variables
  int p; // number of states
  // observations
  vector[r] y[n];
  // system matrices
  // observation equation
  matrix[r, p] F;
  real b;
  // system equation
  matrix[p, p] G;
  real g;
  // initial conditions
  vector[p] m0;
  cov_matrix[p] C0;
}
transformed data {
  matrix[p, p] Ip;
  {
    vector[p] Ip_vector;
    Ip_vector <- rep_vector(1, p);
    Ip <- diag_matrix(Ip_vector);
  }
}
parameters {
  cov_matrix[r] V;  
  vector[2] sigma2;
}
transformed parameters {
  // system variance matrix
  matrix[r, r] W;

  // create observation matrix
  W <- rep_matrix(0, r, r);
  W[2, 2] <- sigma2[1];
  W[3, 3] <- sigma2[2];
}
model {
  y ~ gaussian_dlm_obs(F ', G, V, W, m0, C0);
}
