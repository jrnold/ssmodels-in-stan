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
parameters {
  cov_matrix[p] W;
  cov_matrix[r] V;
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
    matrix[p, p] Ip;
    Ip <- diag_matrix(rep_vector(1, p));

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
      
      // one step ahead predictive distribion of \theta_t | y_{1:(t-1)}
      a[t] <- g + G * m[t];
      R[t] <- quad_form(C[t], G ') + W;
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
      loglik_obs[t] <- - 0.5 * (nobs[t] * log(2 * pi())
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
  // simulated
  vector[p] theta[n + 1];
  // smoothed values
  vector[p] s[n + 1];
  matrix[p, p] S[n + 1];

  // time = n
  // theta_n \sim N(m_n, C_n)
  theta[n + 1] <- multi_normal_rng(m[n + 1], C[n + 1]);
  // time = n - 1, 0
  for (i in 1:n) {
    int t;
    vector[p] h;
    matrix[p, p] H;
    matrix[p, p] Rinv;
    // variable t = time t - 1
    t <- n - i + 1;
    // m[t] = $m_t$ since m is length n + 1
    // C[t] = $C_t$ since C is length n + 1    
    // R[t] = $R_{t + 1}$ since R is length n
    // a[t] = $a_{t + 1}$ since a is length n
    // G[t] = $G_{t + 1}$ since G is length n if time varying
    Rinv <- inverse(R[t]);
    // h_t = m_t + C_t G_{t + 1}' R_{t+1}^{-1} (\theta_{t + 1} - a_{t + 1})
    // For tv G use G[t + 1]
    h <- m[t] + C[t] * G ' * Rinv * (theta[t + 1] - a[t]);
    H <- C[t] - C[t] * G ' * Rinv * G * C[t];
    theta[t] <- multi_normal_rng(h, 0.5 * (H + H '));
  }

  // time = n
  s[n + 1] <- m[n + 1]; 
  S[n + 1] <- C[n + 1];
  // time = n - 1 to 0
  for (i in 1:n) {
    int t;
    matrix[p, p] Rinv;
    matrix[p, p] S_tmp;
    t <- n - i + 1; // i = 1, t = n
    // m[t] = $m_t$ since m is length n + 1
    // C[t] = $C_t$ since C is length n + 1    
    // R[t] = $R_{t + 1}$ since R is length n
    // a[t] = $a_{t + 1}$ since a is length n
    // G[t] = $G_{t + 1}$ since G is length n if time varying
    Rinv <- inverse(R[t]);
    // C_t, G_{t + 1}
    s[t] <- m[t] + C[t] * G ' * Rinv * (s[t + 1] - a[t]);
    S_tmp <- C[t] - C[t] * G ' * Rinv * (R[t] - S[t + 1]) * Rinv * G * C[t];
    S[t] <- 0.5 * (S_tmp ' + S_tmp);
  }
  
}