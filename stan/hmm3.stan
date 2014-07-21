data {
  int n; // number of obs
  int m; // number of states
  vector[n] y;
}
parameters {
  simplex[m] transition[n];
  vector[m] logp[n];
  vector[m] lambda;
}
model {
  // observation
  for (t in 1:n) {
    for (i in 1:m) {
      logp[t, i] <- log_poisson(y[t], lambda[i]);
    }
  }
  // transition
  for (t in 1:n) {
    
  }

  
}
