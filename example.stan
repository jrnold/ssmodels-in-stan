data {
  int N; // number of observations (times)
  int r; // number of variables
  int n; // number of states

  // observations
  vector[r] y[N];
  // observation eq matrices
  matrix[n, r] F;
  vector[r] b;
  // system eq matrices
  vector[n] g;  
  matrix[n, n] G;
  // initial values
  vector[n] m0;
  cov_matrix[n] C0;
}
parameters {
  cov_matrix[r] V;
  cov_matrix[n] W;
}
transformed parameters {
  vector[N] loglik;
  // Kalman filter
  vector[n] m[N + 1];
  matrix[n, n] C[N + 1];
  vector[n] a[N];
  matrix[n, n] C[N];
  vector[r] f;
  matrix[r, r] C[N];

  { // Kalman Filter
    matrix[
    vector[r] err;
    matrix[r, r] Qinv;

    // initial values
    m[1] <- m0;
    C[1] <- C0;
    for (i in 1:N) {
      // prediction equation for theta, p(theta_t | theta_{t-1})
      a[i] <- g + G * m[i];
      R[i] <- G * C[i] * G ' + W;
      // prediction equation for y, p(y_t | y_{1:(t-1)})
      f[i] <- b + F * a[i];
      Q[i] <- F * R[i] * F ' + V;
      // filtered equation p(theta_t | y_{1:t})
      err <- y[i] - f[i];
      Qinv <- inv(Q[i]);
      m[i + 1] <- a[i] + R[i] * F ' Qinv * err;
      C[i + 1] <- R[i] - R[i] * F[i] ' * Qinv * F[i] * R[i]
      // log likelihood
      loglik[i] <- -0.5 * (log(2 * pi()) + log(fabs(Q)) + pow(err, 2) * Qinv);
    }
    increment_log_prob(sum(loglik));
  }

}
model {
}
generated quantities {
  // vector[N + 1] theta;
  // // Backward sample theta
  // theta[N + 1] <- normal_rng(m[N + 1], sqrt(C[N + 1]));
  // for (i in 1:N) {
  //   real h;
  //   real H;
  //   int j;
  //   real m_tmp;
  //   real C_tmp;
  //   j <- N - i;
  //   m_tmp <- m[j + 1];
  //   C_tmp <- C[j + 1];
  //   h <- m_tmp + C_tmp * (theta[j + 2] - a[j + 1]) / R[j + 1];
  //   H <- C_tmp - pow(C_tmp, 2) / R[j + 1];
  //   theta[j + 1] <- normal_rng(h, sqrt(H));
  // }
}
