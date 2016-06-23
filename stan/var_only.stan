functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  int<lower = 1> m;
  int<lower = 1> r;
  vector[1] y[n];
  vector[m] a1;
  cov_matrix[m] P1;
  matrix[m, m] T;
  matrix[1, m] Z;
  matrix[m, r] R;
  vector[m] c;
  vector[1] d;
}
transformed data {
  int filter_sz;
  filter_sz <- ssm_filter_return_size(m, 1);
}
parameters {
  vector<lower = 0.0>[r] sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
  vector[filter_sz] filter_res[n];
  {
    matrix[1, 1] H;
    matrix[r, r] Q;
    Q <- rep_matrix(0.0, r, r);
    for (i in 1:r) {
      Q[i, i] <- pow(sigma_eta[i], 2);
    }
    H <- rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    filter_res <- ssm_filter(y, c, Z, H, d, T, R, Q, a1, P1);
  }
}
model {
  vector[n] ll;
  for (i in 1:n) {
    ll[i] <- ssm_filter_get_loglik(filter_res[i], m, 1);
  }
  print(sum(ll));
  increment_log_prob(sum(ll));
}
