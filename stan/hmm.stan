/* 
   Hidden Markov Model Example

   The following model

   Let $p(\theta_1)$, and for $t = 2, \dots ,T$,
   $$
   \begin{aligned}[t]
   p(\theta_t | \theta_{1:t-1}, \psi) = p(\theta_t | \theta_{t-1}, \psi)  \\
   p(y_t | y_{1:t-1}, \theta_{1:t}) = p(y_t | \theta_t)
   \end{aligned}
   $$

   Stan does not allow sampling discrete parameters. 
   Thus, the sampling procedes in two steps, 

   1. In the ``model`` block, the observed likelihood is calculated by marginalizing out the discrete latent states. 
   2. In the ``generated quantities`` block, the discrete states are sampled conditional on the estimated parameters.  Additionally, statistics of the discrete states are calculated: the conditional probabilities of the states decoding $\Pr(\theta_{t} | y, \psi)$, the global decoding (Viterbi algorithm) $\argmax p(\theta_{1:T} | y_{1:T}, \psi)$.
   
 */ 
data {
  int n; // time
  int m; // number of states
  int y[n];
}
parameters {
  // transition probabilities
  simplex[m] Gamma[m]; // Gamma[from, to]
  // probability parameters
  positive_ordered[m] mu;
  real<lower = 0.0> mu_a;
  real<lower = 0.0> mu_b;  
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
      logp[i, j] <- poisson_log(y[i], mu[j]);
    }
  }
  // calculate the stationary initial distribution
  delta <- (rep_row_vector(1.0, m) /
  	    (diag_matrix(rep_vector(1.0, m))
	    - Gamma_mat + rep_matrix(1.0, m, m))) ';
  
  // HMM Log likelihood
  // Lystig and Hughes (2002) "Exact Computation ..."
  {
    matrix[n, m] log_lambda;
    vector[n] log_Lambda;
    for (j in 1:m) {
      log_lambda[1, j] <- log(delta[j]) + logp[1, j];
    }
    log_Lambda[1] <- log_sum_exp(log_lambda[1]);
    for (t in 2:n) {
      for (j in 1:m) {
	vector[m] tmp;
	for (i in 1:m) {
	  tmp[i] <- (log_lambda[t - 1, i] 
		     + log(Gamma[i, j]));
	}
	log_lambda[t, j] <- (log_sum_exp(tmp) 
			     + logp[t, j] 
			     - log_Lambda[t - 1]);
      }
      log_Lambda[t] <- log_sum_exp(log_lambda[t]);
    }
    llik <- sum(log_Lambda);
  }
}
model {
  increment_log_prob(llik);
  // regularize the states.
  mu ~ gamma(mu_a, mu_b);
}
generated quantities {
  // log forward probabilities
  vector[m] log_alpha[n];
  // log log backward probabilities
  vector[m] log_beta[n];
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
      v <- Gamma_mat * exp(logp[t + 1] + log(phi));
      log_beta[t] <- log(v) + lscale;
      u <- sum(v);
      phi <- v / u;
      lscale <- lscale + log(u);
    }
  }
  // state conditional probabilities
  // uses forward and backward probabilities
  for (i in 1:n) {
    stateprob[i] <- exp(log_alpha[i] + log_beta[i] - llik);
  }
  // Sample values from the states
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
  // Global Decoding (Viterbi algorithm)
  // 
  {
    vector[m] log_xi[n];
    // forwards pass
    log_xi[1] <- log(delta) + logp[1];
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
	log_xi[t, j] <- max_xi_gamma + logp[t, j];
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



