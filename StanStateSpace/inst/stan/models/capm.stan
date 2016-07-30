functions {
  #include ssm.stan
}
data {
  int<lower = 1> n;
  int m;
  vector[m] y[n];
  vector[n] x;

  vector<lower = 0.>[m] a1;
  cov_matrix[m] P1;
  vector<lower = 0.>[m] y_scale;
}
transformed data {
  // system matrices
  vector[m] d[1];
  matrix[m, m] Z[n];
  vector[m] c[1];
  matrix[m, m] T[1];
  matrix[m, m] R[1];
  int filter_sz;
  d[1] = rep_vector(0., m);
  for (t in 1:n) {
    Z[t] = diag_matrix(rep_vector(x[t], m));
  }
  c[1] = rep_vector(0., m);
  T[1] = diag_matrix(rep_vector(1., m));
  R[1] = diag_matrix(rep_vector(1., m));
  # Calculates the size of the vectors returned by ssm_filter
  filter_sz = ssm_filter_size(m, m);
}
parameters {
   cov_matrix[m] Sigma_eta;
   cov_matrix[m] Sigma_epsilon;
/*  vector<lower = 0.>[m] tau_epsilon;
  vector<lower = 0.>[m] tau_eta;
  corr_matrix[m] Rho_epsilon;
  corr_matrix[m] Rho_eta;*/
}
transformed parameters {
  // cov_matrix[m] Sigma_eta;
  // cov_matrix[m] Sigma_epsilon;
  // Sigma_epsilon = quad_form_diag(Rho_epsilon, tau_epsilon);
  // Sigma_eta = quad_form_diag(Rho_eta, tau_eta .* tau_epsilon);
}
model {
  matrix[m, m] H[1];
  matrix[m, m] Q[1];
  H = rep_array(Sigma_epsilon, 1);
  Q = rep_array(Sigma_eta, 1);
  y ~ ssm_lpdf(d, Z, H,
               c, T, R, Q,
               a1, P1);
  // Rho_epsilon ~ lkj_corr(5.);
  // Rho_eta ~ lkj_corr(5.);
  // tau_epsilon ~ cauchy(0., y_scale);
  // tau_eta ~ cauchy(0., 1.);
}
generated quantities {
  vector[filter_sz] filtered[n];
  vector[m] beta[n];
  vector[m] beta_mean[n];
  {
    matrix[m, m] H[1];
    matrix[m, m] Q[1];
    H = rep_array(Sigma_epsilon, 1);
    Q = rep_array(Sigma_eta, 1);
    // filtering
    filtered = ssm_filter(y, d, Z, H,
                          c, T, R, Q,
                          a1, P1);
    // mean states
    beta_mean = ssm_smooth_states_mean(filtered, Z, c, T, R, Q);
    // sampling states
    beta = ssm_simsmo_states_rng(filtered,
                                 d, Z, H,
                                 c, T, R, Q,
                                 a1, P1);

  }
}
