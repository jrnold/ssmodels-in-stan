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
  cov_matrix[p] W;
  vector<lower = 0.0>[r] V;
}
transformed parameters {
  // log-likelihood
  vector[r] loglik_obs[n];
  real loglik;
  // prior of state: p(theta_t | y_t, ..., y_{t-1})
  vector[p] a[n];
  cov_matrix[p] R[n];
  // likelihood of obs: p(y_t | y_t, ..., y_t-1)
  vector[r] f[n];
  vector<lower=0.0>[r] Q[n];
  // posterior of states: p(theta_t | y_t, ..., y_t)
  vector[p] m[n + 1];
  cov_matrix[p] C[n + 1];
  // create observation matrix
  {
    real err;
    vector[p] K;
    matrix[p, p] J;
    vector[p] m_tmp;
    matrix[p, p] R_tmp;
    matrix[p, p] C_tmp;
    vector[p] Fj;
    // set initial states
    m[1] <- m0;
    C[1] <- C0;
    // loop through observations
    for (t in 1:n) {
      a[t] <- g + G * m[t];
      R_tmp <- quad_form(C[t], G ') + W;
      R[t] <- 0.5 * (R_tmp + R_tmp ');
      m_tmp <- a[t];
      C_tmp <- R[t];
      ## filter by observation
      for (j in 1:r) {
	Fj <- row(F, j) ';
	// one step ahead predictive distribion of \theta_t | y_{1:(t-1)}
	// one step ahead predictive distribion of y_t | y_{1:(t-1)}
	f[t, j] <- b[j] + dot_product(Fj, m_tmp);
	Q[t, j] <- Fj ' * C_tmp * Fj + V[j]; 
	// forecast error
	err <- y[t, j] - f[t, j];
	// Kalman gain
	K <- C_tmp * Fj / Q[t, j];
	// posterior distribution of \theta_t | y_{1:t}
	m_tmp <- m_tmp + K * err;
	// matrix used in Joseph stabilized form
	J <- (Ip - K * Fj ');
	C_tmp <- quad_form(C_tmp, J ') + K ' * K * V[j];
	// log likelihood
	loglik_obs[t, j] <- - 0.5 * (log(2 * pi())
				     + log(Q[t, j])
				     + pow(err, 2) / Q[t, j]);
      }
      m[t + 1] <- m_tmp;
      C[t + 1] <- 0.5 * (C_tmp + C_tmp ');
    }
  }
  loglik <- 0.0;
  for (i in 1:n) {
    loglik <- loglik + sum(loglik_obs[i]);
  }
}
model {
  increment_log_prob(loglik);
}
generated quantities {
  vector[p] theta[n + 1];
  theta[n + 1] <- multi_normal_rng(m[n + 1], C[n + 1]);
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