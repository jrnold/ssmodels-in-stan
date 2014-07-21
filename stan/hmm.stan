data {
  int n; // time
  int m; // number of states
  int y[n];
  // vector[m] lambda;
}
parameters {
  // transition probabilities
  simplex[m] Gamma[m]; // Gamma[from, to]
  // probability parameters
  positive_ordered[m] lambda;
  real<lower = 0.0> lambda_a;
  real<lower = 0.0> lambda_b;  
}
transformed parameters {
  vector<upper = 0>[m] logp[n];
  real llik;
  
  simplex[m] delta;
  matrix[m, m] Gamma_mat;
  for (i in 1:m) {
    for (j in 1:m) {
      Gamma_mat[i, j] <- Gamma[i, j];
    }
  }
  
  for (i in 1:n) {
    for (j in 1:m) {
      logp[i, j] <- poisson_log(y[i], lambda[j]);
    }
  }
  // calculate the stationary initial distribution
  delta <- (rep_row_vector(1.0, m) /
  	    (diag_matrix(rep_vector(1.0, m))
	    - Gamma_mat + rep_matrix(1.0, m, m))) ';
  
  // HMM Log likelihood
  // Zucchini and MacDonald, p. 47,
  {
    vector[m] phi;
    real u;
    vector[m] v;
    llik <- 0;
    for (i in 1:n) {
      if (i == 1) {
      	v <- exp(log(delta) + logp[1]);
      } else {
	v <- exp(log(phi ' * Gamma_mat) + logp[i] ') ';
      }
      u <- sum(v);
      phi <- v / u;
      llik <- llik + log(u);
    }
  }
}
model {
  increment_log_prob(llik);
  lambda ~ gamma(lambda_a, lambda_b);
}
generated quantities {
  vector[m] log_alpha[n];  
  vector[m] log_beta[n];
  vector[m] stateprob[n];
  int states[n];

  // forward probailities
  {
    real u;
    vector[m] phi;
    vector[m] v;
    real lscale;
    phi <- delta;
    lscale <- 0.0;
    for (i in 1:n) {
      if (i == 1) {
      	v <- exp(log(delta) + logp[1]);
      } else {
  	v <- exp(log(phi ' * Gamma_mat) + logp[i] ') ';
      }
      u <- sum(v);
      lscale <- lscale + log(u);
      phi <- v / u;
      log_alpha[i] <- log(phi) + lscale;
    }
  }
  // backward probilities
  {
    int t;
    real u;
    vector[m] phi;
    vector[m] v;
    real lscale;

    log_beta[n] <- rep_vector(0.0, m);
    phi <- rep_vector(1.0 / m, m);
    lscale <- log(m);
    for (i in 1:(n - 1)) {
      t <- n - i;
      v <- Gamma_mat * exp(logp[t + 1] + log(phi)) ;
      log_beta[t] <- log(v) + lscale;
      u <- sum(v);
      phi <- v / u;
      lscale <- lscale + log(u);
    }
  }
  // state conditional probabilities
  for (i in 1:n) {
    stateprob[i] <- exp(log_alpha[i] + log_beta[i] - llik);
  }
  // sampling states
  states[n] <- categorical_rng(exp(log_alpha[n]) / sum(exp(log_alpha[n])));
  for (i in 1:(n - 1)) {
    int t;
    vector[m] theta;
    t <- n - i;
    for (j in 1:m) {
      theta[j] <- exp(log_alpha[t, j]) * Gamma[j, states[t + 1]];
    }
    theta <- theta / sum(theta);
    states[t] <- categorical_rng(theta);
  }
  // viterbi
  // TODO
}



