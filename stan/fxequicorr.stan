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
  vector[r] b;
  // system equation
  matrix[p, p] G;
  vector[p] g;
  // initial conditions
  vector[p] m0;
  cov_matrix[p] C0;

}
transformed data {
  matrix[p, p] Ip;
  Ip <- diag_matrix(rep_vector(1, p));
}
parameters {
  vector<lower = 0.0>[p] Wdiag;
  vector<lower = 0.0>[r] Vdiag;
  real<lower = -1.0, upper = 1.0> rho;
}
transformed parameters {
  // system variance matrix
  cov_matrix[p] W;
  cov_matrix[r] V;

  // log-likelihood
  real loglik_obs[n];
  real loglik;
  // prior of state: p(theta_t | y_t, ..., y_{t-1})
  vector[p] a[n];
  matrix[p, p] R[n];
  // likelihood of obs: p(y_t | y_t, ..., y_t-1)
  vector[r]  f[n];
  matrix[r, r] Q[n];
  // posterior of states: p(theta_t | y_t, ..., y_t)
  vector[p] m[n + 1];
  matrix[p, p] C[n + 1];

  W <- diag_matrix(Wdiag);
  for (i in 1:r) {
    for (j in 1:r) {
      if (i == j) {
	V[i, i] <- pow(Vdiag[i], 2);
      } else {
	V[i, j] <- Vdiag[i] * Vdiag[j] * rho;
      }
    }
  }
  {

    // set initial states
    m[1] <- m0;
    C[1] <- C0;
    
    // loop through observations    
    for (t in 1:n) {
      vector[r] err;
      matrix[p, r] K;
      matrix[r, r] Q_tmp;      
      matrix[r, r] Qinv;
      matrix[p, p] J;
      matrix[p, p] R_tmp;
      
      // one step ahead predictive distribion of \theta_t | y_{1:(t-1)}
      a[t] <- g + G * m[t];
      R_tmp <- quad_form(C[t], G ') + W;
      R[t] <- 0.5 * (R_tmp + R_tmp ');
      // one step ahead predictive distribution of y_t | y_{1:(t-1)}
      f[t] <- b + F * a[t];
      Q_tmp <- quad_form(R[t], F ') + V;
      Q[t] <- 0.5 * (Q_tmp + Q_tmp ');
      // forecast error
      err <- y[t] - f[t];
      // Kalman gain
      Qinv <- inverse_spd(Q[t]);
      K <- R[t] * F ' * Qinv;
      // posterior distribution of \theta_t | y_{1:t}
      m[t + 1] <- a[t] + K * err;
      // matrix used in Joseph stabilized form
      J <- (Ip - K * F);
      C[t + 1] <- quad_form(R[t], J ') + quad_form(V, K ');
      // log likelihood
      loglik_obs[t] <- - 0.5 * (r * log(2 * pi())
				+ log_determinant(Q_tmp)
				+ quad_form(Qinv, err));
    }
  }
  loglik <- sum(loglik_obs);
}
model {
  increment_log_prob(loglik);
}
generated quantities {
  vector[p] theta[n + 1];
  theta[n + 1] <- multi_normal_rng(m[n + 1], C[n + 1]);
  // iterate backwards over observations
  for (i in 1:n) {
    int t;
    vector[p] h;
    matrix[p, p] H;
    matrix[p, p] Rinv;
    t <- n - i + 1;
    Rinv <- inverse(R[t]);
    // sample 
    h <- m[t] + C[t] * G ' * Rinv * (theta[t + 1] - a[t]);
    H <- C[t] - C[t] * G ' * Rinv * G * C[t];
    theta[t] <- multi_normal_rng(h, 0.5 * (H + H '));
  }
}