data {
  int T; // time
  int s; // number of states
  vector y[T];
}
parameters {
  // transition probabilities
  // TODO: switch to array of unitary vectors
  simplex[s] gamma[s]; // [from, to]
  // initial state probabilities
  simplex[s] delta;
  // log observation probabilities
  matrix<upper = 0>[T, s] log_pi;
}
transformed parameters {
  matrix[T, s] log_lambda;
  vector[T] log_Lambda;
  for (j in 1:s) {
    log_lambda[1, j] <- log(delta[j]) + log_pi[1, j];
  }
  log_Lambda[1] <- log_sum_exp(lambda[1]);
  for (t in 2:T) {
    for (j in 1:s) {
      vector[s] tmp;
      for (i in 1:s) {
	tmp[i] <- (log_lambda[t - 1, i] 
		   + log(gamma[i, j]));
      }
      log_lambda[j] <- (log_sum_exp(tmp) 
			+ log_pi[t, j] 
			- log_Lambda[t - 1]);
    }
    log_Lambda[t] <- log_sum_exp(log_lambda[t]);
  }
}
model {
  increment_log_prob(sum(log_Lambda));
}
generated quantities {
  // backward algorithm
}



