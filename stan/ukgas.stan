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
  {
    vector[p] Ip_vector;
    Ip_vector <- rep_vector(1, p);
    Ip <- diag_matrix(Ip_vector);
  }
}
parameters {
  cov_matrix[r] V;  
  vector<lower = 0.0>[2] sigma2;
}
transformed parameters {
  // system variance matrix
  matrix[p, p] W;
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
  // create observation matrix
  W <- rep_matrix(0, p, p);
  W[2, 2] <- sigma2[1];
  W[3, 3] <- sigma2[2];
  {
    vector[r] err;
    matrix[p, r] K;
    matrix[r, r] Qinv;
    matrix[p, p] J;
    // set initial states
    m[1] <- m0;
    C[1] <- C0;
    // loop through observations
    for (t in 1:n) {
      // one step ahead predictive distribion of \theta_t | y_{1:(t-1)}
      a[t] <- g + G * m[t];
      R[t] <- quad_form(C[t], G ') + W;
      // one step ahead predictive distribion of y_t | y_{1:(t-1)}
      f[t] <- b + F * a[t];
      Q[t] <- quad_form(R[t], F ') + V;      
      // forecast error
      err <- y[t] - f[t];
      // Kalman gain
      Qinv <- inverse(Q[t]);
      K <- R[t] * F ' * Qinv;
      // posterior distribution of \theta_t | y_{1:t}
      m[t + 1] <- a[t] + K * err;
      // matrix used in Joseph stabilized form
      J <- (Ip - K * F);
      C[t + 1] <- quad_form(R[t], J ') + quad_form(V, K ');
      // log likelihood
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
generated quantities {
  vector[p] theta[n + 1];
  vector[p] s[n + 1];
  matrix[p, p] S[n + 1];  
  theta[n + 1] <- multi_normal_rng(m[n + 1], C[n + 1]);
  s[n + 1] <- m[n + 1];
  S[n + 1] <- C[n + 1];
  for (i in 1:n) {
    int t;
    vector[p] h;
    matrix[p, p] H;
    matrix[p, p] Rinv;
    matrix[p, p] S_tmp;
    t <- n - i + 1;
    Rinv <- inverse(R[t]);
    // smoother
    s[t] <- m[t] + C[t] * G ' * Rinv * (s[t + 1] - a[t]);
    S_tmp <- C[t] - C[t] * G ' * Rinv * (R[t] - S[t + 1]) * Rinv * G * C[t];
    S[t] <- 0.5 * (S_tmp + S_tmp ');
    // sample 
    h <- m[t] + C[t] * G ' * Rinv * (theta[t + 1] - a[t]);
    H <- C[t] - C[t] * G ' * Rinv * G * C[t];
    theta[t] <- multi_normal_rng(h, 0.5 * (H + H '));
  }
}