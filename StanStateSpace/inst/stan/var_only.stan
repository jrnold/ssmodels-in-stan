functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  int<lower = 1> m;
  int<lower = 1> q;
  vector[1] y[n];
  vector[m] a1;
  cov_matrix[m] P1;
  matrix[m, m] T;
  matrix[1, m] Z;
  matrix[m, q] R;
  vector[m] c;
  vector[1] d;
}
transformed data {
  int filter_sz;
  int p;
  p = 1;
  filter_sz = ssm_filter_size(m, p);
}
parameters {
  vector<lower = 0.0>[q] sigma_eta;
  real<lower = 0.0> sigma_epsilon;
}
transformed parameters {
}
model {
  {
    matrix[1, 1] H;
    matrix[q, q] Q;
    Q = rep_matrix(0.0, q, q);
    for (i in 1:q) {
      Q[i, i] = pow(sigma_eta[i], 2);
    }
    H = rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    target += ssm_constant_lpdf(y | c, Z, H, d, T, R, Q, a1, P1);
  }
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[q + q * q] eta[n];
  vector[p + p * p] eps[n];
  vector[m + m * m] alpha[n];
  vector[m] alpha2[n];
  vector[m + m * m] alpha3[n];
  vector[2 * p + m + q] sims[1];
  {
    matrix[1, 1] H;
    matrix[q, q] Q;
    Q = rep_matrix(0.0, q, q);
    for (i in 1:q) {
      Q[i, i] = pow(sigma_eta[i], 2);
    }
    H = rep_matrix(pow(sigma_epsilon, 2), 1, 1);
    filtered = ssm_filter(y, rep_array(d, 1), rep_array(Z, 1), rep_array(H, 1),
                             rep_array(c, 1), rep_array(T, 1), rep_array(R, 1),
                             rep_array(Q, 1), a1, P1);
    eta = ssm_smooth_eta(filtered,
                         rep_array(Z, 1), rep_array(T, 1), rep_array(R, 1),
                         rep_array(Q, 1));
    eps = ssm_smooth_eps(filtered, rep_array(Z, 1), rep_array(H, 1), rep_array(T, 1));
    alpha = ssm_smooth_state(filtered, rep_array(Z, 1), rep_array(T, 1));
    alpha2 = ssm_smooth_faststate(filtered, rep_array(c, 1),
                                  rep_array(Z, 1), rep_array(T, 1), rep_array(R, 1),
                                  rep_array(Q, 1));
    alpha3 = ssm_filter_states(filtered, rep_array(Z, 1));
    sims = ssm_sim_rng(1, rep_array(c, 1), rep_array(Z, 1), rep_array(H, 1),
                       rep_array(d, 1), rep_array(T, 1), rep_array(R, 1),
                       rep_array(Q, 1), a1, P1);

  }
}
