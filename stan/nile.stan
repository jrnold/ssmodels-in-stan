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
  cov_matrix[p] W;
}
transformed parameters {
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

  {
    vector[r] err;
    matrix[p, r] K;
    matrix[r, r] Qinv;
    matrix[p, p] C_tmp;
    matrix[p, p] J;
    
    // set initial states
    m[1] <- m0;
    C[1] <- C0;
    // loop through observations
    for (t in 1:n) {
      // one step ahead predictive distribion of \theta_t | y_{1:(t-1)}
      a[t] <- g + G * m[t];
      // R[t] <- G * C[t] * G ' + W;
      R[t] <- quad_form(C[t], G ') + W;
      // one step ahead predictive distribion of y_t | y_{1:(t-1)}
      f[t] <- b + F * a[t];
      // Q[t] <- F * R[t] * F ' + V;
      Q[t] <- quad_form(R[t], F ') + V;      
      Qinv <- inverse(Q[t]);
      // error
      err <- y[t] - f[t];
      // Kalman gain
      K <- R[t] * F ' * Qinv;
      // posterior distribution of \theta_t | y_{1:t}
      m[t + 1] <- a[t] + K * err;
      // matrix used in Joseph stabilized form
      // C_tmp <- R[t] - K * Q[t] * K ';
      // C_tmp <- R[t] - quad_form(Q[t], K ');
      // C[t + 1] <- 0.5 * (C_tmp + C_tmp ');
      J <- (Ip - K * F);
      C[t + 1] <- quad_form(R[t], J ') + quad_form(V, K ');
      // log likelihood
      // loglik_obs[t] <- - 0.5 * (r * log(2 * pi())
      // 				+ log_determinant(Q[t])
      // 				+ err ' * Qinv * err);
      loglik_obs[t] <- - 0.5 * (r * log(2 * pi())
      				+ log_determinant(Q[t])
      				+ quad_form(Qinv, err));

    }
  }
  loglik <- sum(loglik_obs);
}
model {
  increment_log_prob(loglik);
}
