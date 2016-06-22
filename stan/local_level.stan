functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  vector[1] y[n];
  vector<lower = 0.0>[1] a1;
  cov_matrix[1] P1;
}
transformed data {
  matrix[1, 1] T;
  matrix[1, 1] Z;
  matrix[1, 1] R;
  vector[1] c;
  vector[1] d;
  T[1, 1] <- 1.0;
  Z[1, 1] <- 1.0;
  R[1, 1] <- 1.0;
  c[1] <- 0.0;
  d[1] <- 0.0;
}
parameters {
  real<lower = 0.0> sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
  vector[6] filter_res[n];
  {
    matrix[1, 1] H;
    matrix[1, 1] Q;
    H <- rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    Q <- rep_matrix(pow(sigma_eta, 2), 1, 1);
    filter_res <- ssm_filter(y, c, Z, H, d, T, R, Q, a1, P1);
  }
}
model {
  vector[n] ll;
  for (i in 1:n) {
    ll[i] <- ssm_filter_get_loglik(filter_res[i], 1, 1);
  }
  increment_log_prob(sum(ll));
}
