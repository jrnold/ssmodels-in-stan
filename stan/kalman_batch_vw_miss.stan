data {
  // dimensions
  int n; // number of observations
  int r; // number of variables
  int p; // number of states
  // observations
  vector[r] y[n];
  // number of observed variables per observation
  int<lower = 0, upper = r> y_obs_n[n];
  // indices of the observed variables per observation
  // for row t, the code will only use the t entries
  int<lower = 0, upper = r> y_obs_i[n, r];
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
  real<lower=0.0> Qmiss;
}
parameters {
  cov_matrix[p] W;
  cov_matrix[r] V;
}
transformed parameters {
  // log-likelihood by observation
  real loglik_obs[n];
  // prior state: p(theta_t | y_t, ..., y_{t-1})
  vector[p] a[n];
  matrix[p, p] R[n];
  // marginal likelihood: p(y_t | y_t, ..., y_t-1)
  vector[r] f[n];
  matrix[r, r] Q[n];
  // state posterior: p(theta_t | y_t, ..., y_t)
  vector[p] m[n + 1];
  matrix[p, p] C[n + 1];

  // BEGIN KALMAN FILTER
  {
    matrix[p, p] Ip;
    Ip <- diag_matrix(rep_vector(1, p));

    // set initial states
    m[1] <- m0;
    C[1] <- C0;
    
    // loop through observations    
    for (t in 1:n) {
      vector[nobs[t]] y_tmp;
      matrix[nobs[t], p] F_tmp;
      vector[nobs[t]] b_tmp;
      vector[nobs[t], nobs[t]] V_tmp;
      vector[nobs[t]] f_tmp;      
      matrix[nobs[t], nobs[t]] Q_tmp;

      vector[nobs[t]] err;
      matrix[p, nobs[t]] K;
      matrix[nobs[t], nobs[t]] Qinv;
      matrix[p, p] J;
      
      // one step ahead predictive distribion of \theta_t | y_{1:(t-1)}
      a[t] <- g + G * m[t];
      R[t] <- quad_form(C[t], G ') + W;
      f[t] <- rep_vector(0.0, r);
      Q[t] <- diag_matrix(rep_vector(Qmiss, r));
      if (nobs[t] == 0) {
	// if all observations missing
	m[t + 1] <- a[t];
	R[t + 1] <- R[t];
	loglik_obs[t] <- 0.0;
      } else {
	if (nobs[t] == r) {
	  // if no observations missing
	  y_tmp <- y;
	  b_tmp <- b;
	  F_tmp <- F;
	  V_tmp <- V;
	} else {
	  // if at least one observation missing
	  for (i in 1:nobs[t]) {
	    y_tmp[i] <- y[obs_ind[i]];
	    b_tmp[i] <- y[obs_ind[i]];
	    for (j in 1:p) {
	      F[i, j] <- F[obs_ind[i], j];
	    }
	    for (j in 1:nobs[t]) {
	      V_tmp[i, j] <- V[obs_ind[i], obs_ind[j]];
	    }
	  }
	}
	// one step ahead predictive distribution of y_t | y_{1:(t-1)}
	f_tmp <- b_tmp + F_tmp * a[t];
	Q_tmp <- quad_form(R[t], F_tmp ') + V_tmp;
	// forecast error
	err <- y_tmp[t] - f_tmp[t];
	// Kalman gain
	Qinv <- inverse_spd(Q[t]);
	K <- R[t] * F_tmp ' * Qinv;
	// posterior distribution of \theta_t | y_{1:t}
	// calculated using Joseph stabilized form
	m[t + 1] <- a[t] + K * err;
	J <- (Ip - K * F_tmp);
	C[t + 1] <- quad_form(R[t], J ') + quad_form(V_tmp, K ');
	// log likelihood of the observation
	loglik_obs[t] <- - 0.5 * (nobs[t] * log(2 * pi())
				  + log_determinant(Q_tmp)
				  + quad_form(Qinv, err));
	if (nobs[t] == r) {
	  f[t] <- f_tmp;
	  Q[t] <- Q_tmp;
	} else {
	  for (i in 1:nobs[t]) {
	    f[obs_ind[i]] <- f_tmp[i];
	    for (j in 1:nobs[t]) {
	      Q[obs_ind[i], obs_ind[j]] <- Q_tmp[i, j];
	    }
	  }
	}
      }
    }
  }
  // END KALMAN FILTER
}
model {
  increment_log_prob(sum(loglik_obs));
}
generated quantities {
  // sample of theta 
  vector[p] theta[n + 1];
  // BEGIN BACKWARD SAMPLING
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
  // END BACKWARD SAMPLING
}