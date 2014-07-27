/* 
   Chib (1998) Changepoint Model for Coal Mining Disaster Data
   
 */ 
data {
  int n; // time
  int m; // number of changepoints
  int y[n];
  // priors
  real<lower = 0.0> rho_a;
  real<lower = 0.0> rho_b;
  real<lower = 0.0> mu_a;
  real<lower = 0.0> mu_b;
}
transformed data {
  simplex[m] delta;  
  delta <- rep_vector(0.0, m);
  delta[1] <- 1.0;
}
parameters {
  // transition probabilities
  vector<lower = 0.0, upper = 1.0>[m - 1] rho;
  // observation probability parameters
  vector<lower = 0.0>[m] mu;
}
transformed parameters {
  matrix[m, m] Gamma; // transition matrix
  vector<lower = 0, upper=1>[m] p[n];
  real llik;

  Gamma <- rep_matrix(0.0, m, m);
  for (i in 1:(m - 1)) {
    Gamma[i, i] <- rho[i];
    Gamma[i, i + 1] <- 1.0 - rho[i];
  }
  Gamma[m, m] <- 1.0;

  for (i in 1:n) {
    for (j in 1:m) {
      p[i, j] <- exp(poisson_log(y[i], mu[j]));
    }
  }
  // HMM Log likelihood
  // Lystig and Hughes (2002) "Exact Computation ..."
  {
    matrix[n, m] lambda;
    vector[n] Lambda;
    for (j in 1:m) {
      // cannot use logs because delta may be 0.
      lambda[1, j] <- delta[j] * p[1, j];
    }
    Lambda[1] <- sum(lambda[1]);
    for (t in 2:n) {
      for (j in 1:m) {
	vector[m] tmp;
	for (i in 1:m) {
	  tmp[i] <- lambda[t - 1, i] * Gamma[i, j];
	}
	lambda[t, j] <- (sum(tmp) 
			 * p[t, j] 
			 / Lambda[t - 1]);
      }
      Lambda[t] <- sum(lambda[t]);
    }
    llik <- sum(log(Lambda));
  }
}
model {
  increment_log_prob(llik);
  mu ~ gamma(mu_a, mu_b);
  rho ~ beta(rho_a, rho_b);
}
generated quantities {
  // log forward probabilities
  vector<upper=0>[m] log_fwd_prob[n];
  // log log backward probabilities
  vector<upper=0>[m] log_bckwd_prob[n];
  // state probabilities
  vector[m] stateprob[n];
  // global decoding, calcualated with Viterbi algorithm
  int viterbi[n];
  // sample of latent states
  int states[n];

  // Forward probabilities
  // this could be made more efficient by reusing forward probabilities
  // calculated in the likelihood
  {
    real u;
    vector[m] phi;
    vector[m] v;
    real lscale;
    phi <- delta;
    lscale <- 0.0;
    for (i in 1:n) {
      if (i == 1) {
      	v <- delta .* p[1] ;
      } else {
  	v <- (phi ' * Gamma .* p[i] ') ';
      }
      u <- sum(v);
      lscale <- lscale + log(u);
      phi <- v / u;
      log_fwd_prob[i] <- log(phi) + lscale;
    }
  }
  // backward probilities
  {
    int t;
    real u;
    vector[m] phi;
    vector[m] v;
    real lscale;

    log_bckwd_prob[n] <- rep_vector(0.0, m);
    phi <- rep_vector(1.0 / m, m);
    lscale <- log(m);
    for (i in 1:(n - 1)) {
      t <- n - i;
      v <- Gamma * diag_matrix(p[t + 1]) * phi;
      log_bckwd_prob[t] <- log(v) + lscale;
      u <- sum(v);
      phi <- v / u;
      lscale <- lscale + log(u);
    }
  }
  // state conditional probabilities
  // uses forward and backward probabilities
  for (i in 1:n) {
    stateprob[i] <- exp(log_fwd_prob[i] + log_bckwd_prob[i] - llik);
  }
  // Sample values from the states
  states[n] <- categorical_rng(exp(log_fwd_prob[n] - log_sum_exp(log_fwd_prob[n])));
  for (i in 1:(n - 1)) {
    int t;
    vector[m] theta;
    t <- n - i;
    for (j in 1:m) {
      theta[j] <- exp(log_fwd_prob[t, j]) * Gamma[j, states[t + 1]];
    }
    theta <- theta / sum(theta);
    states[t] <- categorical_rng(theta);
  }
  // Global Decoding (Viterbi algorithm)
  // 
  {
    vector[m] log_xi[n];
    // forwards pass
    log_xi[1] <- log(delta) + log(p[1]);
    for (t in 2:n) {
      for (j in 1:m) {
  	real max_xi_gamma;
  	real tmp;
  	max_xi_gamma <- negative_infinity();
  	for (i in 1:m) {
  	  tmp <- log_xi[t - 1, i] + log(Gamma[i, j]);
  	  if (tmp > max_xi_gamma) {
  	    max_xi_gamma <- tmp;
  	  }
  	}
  	log_xi[t, j] <- max_xi_gamma + log(p[t, j]);
      }
    }
    // backwards pass
    {
      real tmp;
      tmp <- negative_infinity();
      for (i in 1:m) {
    	if (log_xi[n, i] > tmp) {
  	  tmp <- log_xi[n, i];
    	  viterbi[n] <- i;
    	}
      }
      for (s in 1:(n - 1)) {
    	int t;
    	real tmp1;
    	t <- n - s;
    	tmp1 <- negative_infinity();
    	for (i in 1:m) {
    	  real tmp2;	
    	  tmp2 <- log_xi[t, i] + log(Gamma[i, viterbi[t + 1]]);
    	  if (tmp2 > tmp1) {
  	    tmp1 <- tmp2;
    	    viterbi[t] <- i;
    	  }
    	}
      }
    }
  }
}

