/////////// Utility Functions  /////////////////
matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}

matrix to_matrix_colwise(vector v, int m, int n) {
  matrix[m, n] res;
  for (j in 1:n) {
    for (i in 1:m) {
      res[i, j] <- v[(j - 1) * m + m];
    }
  }
  return res;
}

matrix to_matrix_rowwise(vector v, int m, int n) {
  matrix[m, n] res;
  for (i in 1:n) {
    for (j in 1:m) {
      res[i, j] <- v[(i - 1) * n + n];
    }
  }
  return res;
}

vector to_vector_colwise(matrix x) {
  vector[num_elements(x)] res;
  int n;
  int m;
  n <- rows(x);
  m <- cols(x);
  for (i in 1:n) {
    for (j in 1:m) {
      res[n * (j - 1) + i] <- x[i, j];
    }
  }
  return res;
}

vector to_vector_rowwise(matrix x) {
  vector[num_elements(x)] res;
  int n;
  int m;
  n <- rows(x);
  m <- cols(x);
  for (i in 1:rows(x)) {
    for (j in 1:cols(x)) {
      res[(i - 1) * m + j] <- x[i, j];
    }
  }
  return res;
}

/////////// SSM Utilities ////////////////



/////////// SSM Filter  /////////////////

// Length of vectors that SSM returns
// value    size      location
// log-lik  1         1
// v        p         2
// F^-1     p^2       2 + p
// K        mp        2 + p + p^2
// a_t      m         2 + p + p^2 + mp
// P^t      m * m     2 + p + p^2 + mp + m

// rows (loglik, v, Finv, K, a, P)
// cols (loc, length)
int[,] ssm_filter_idx(int m, int p) {
  int sz[6, 3];
  // loglike
  sz[1, 1] <- 1;
  // v
  sz[2, 1] <- p;
  // Finv
  sz[3, 1] <- p * p;
  // K
  sz[4, 1] <- m * p;
  // a
  sz[5, 1] <- m;
  // P
  sz[6, 1] <- m * m;

  // Fill in start and stop points
  sz[1, 2] <- 1;
  sz[1, 3] <- sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:6) {
    sz[i, 2] <- sz[i - 1, 3] + 1;
    sz[i, 3] <- sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}

int ssm_filter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx <- ssm_filter_idx(m, p);
  sz <- idx[6, 3];
  return sz;
}

real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y <- x[1];
  return y;
}

vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] y;
  y <- segment(x, 2, p);
  return y;
}

matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p, p] y;
  y <- to_matrix_colwise(segment(x, 2 + p, p * p), p, p);
  return y;
}

matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  y <- to_matrix_colwise(segment(x, 2 + p + p * p, m * p), m, p);
  return y;
}

vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  y <- segment(x, 2 + p + p * p + m * p, m);
  return y;
}

matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  y <- to_matrix_colwise(segment(x,  2 + p + p * p + m * p + m, m * m), m, m);
  return y;
}


// Filtering and Log-likelihood
vector ssm_filter_update_a(vector a, vector c, matrix T, vector v, matrix K) {
  vector[num_elements(a)] a_new;
  a_new <- T * a + K * v + c;
  return a_new;
}

matrix ssm_filter_update_P(matrix P, matrix Z, matrix T,
                           matrix Q, matrix R, matrix K) {
  matrix[rows(P), cols(P)] P_new;
  P_new <- to_symmetric_matrix(T * P * (T - K * Z)' + quad_form(Q, R'));
  return P_new;
}

vector ssm_filter_update_v(vector y, vector a, vector d, matrix Z) {
  vector[num_elements(y)] v;
  v <- y - Z * a - d;
  return v;
}

matrix ssm_filter_update_F(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] F;
  F <- quad_form(P, Z') + H;
  return F;
}

matrix ssm_filter_update_Finv(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] Finv;
  Finv <- inverse(ssm_filter_update_F(P, Z, H));
  return Finv;
}

matrix ssm_filter_update_K(matrix P, matrix T, matrix Z, matrix Finv) {
  matrix[cols(Z), rows(Z)] K;
  K <- T * P * Z' * Finv;
  return K;
}

matrix ssm_filter_update_L(matrix Z, matrix T, matrix K) {
  matrix[rows(T), cols(T)] L;
  L <- T - K * Z;
  return L;
}

real ssm_filter_update_ll(vector v, matrix Finv) {
  real ll;
  int p;
  p <- num_elements(v);
  // det(A^{-1}) = 1 / det(A) -> log det(A^{-1}) = - log det(A)
  ll <- (- 0.5 * (p * log(2 * pi())
         - log_determinant(Finv)
         + quad_form(Finv, v)));
  return ll;
}

real ssm_lpdf(vector[] y,
                vector c, matrix Z, matrix H,
                vector d, matrix T, matrix R, matrix Q,
                vector a1, matrix P1) {
  real ll;
  int n;
  int m;
  int p;

  n <- size(y); // number of obs
  m <- cols(Z);
  p <- rows(Z);
  {
    vector[n] ll_obs;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;

    a <- a1;
    P <- P1;
    for (t in 1:n) {
      v <- ssm_filter_update_v(y[t], a, d, Z);
      Finv <- ssm_filter_update_Finv(P, Z, H);
      K <- ssm_filter_update_K(P, T, Z, Finv);
      ll_obs[t] <- ssm_filter_update_ll(v, Finv);
      // don't save a, P for last iteration
      if (t < n) {
        a <- ssm_filter_update_a(a, c, T, v, K);
        P <- ssm_filter_update_P(P, Z, T, Q, R, K);
      }
    }
    ll <- sum(ll_obs);
  }
  return ll;
}

void ssm_lp(vector[] y,
            vector c, matrix Z, matrix H,
            vector d, matrix T, matrix R, matrix Q,
            vector a1, matrix P1) {
  real ll;
  ll <- ssm_lpdf(y, c, Z, H, d, T, R, Q, a1, P1);
  increment_log_prob(ll);
}

// Filtering
vector[] ssm_filter(vector[] y,
                    vector c, matrix Z, matrix H,
                    vector d, matrix T, matrix R, matrix Q,
                    vector a1, matrix P1) {

  // returned data
  vector[ssm_filter_size(cols(Z), rows(Z))] res[size(y)];
  int q;
  int n;
  int p;
  int m;

  // sizes
  n <- size(y); // number of obs
  p <- rows(Z); // obs size
  m <- cols(Z); // number of states
  q <- cols(Q); // number of state disturbances

  //print("Sizes: n = ", m, ", p = ", n, ", m = ", m, ", q = ", q);
  {
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    real ll;
    int idx[6, 3];

    idx <- ssm_filter_idx(m, p);
    a <- a1;
    P <- P1;
    for (t in 1:n) {
      // updating
      v <- ssm_filter_update_v(y[t], a, d, Z);
      Finv <- ssm_filter_update_Finv(P, Z, H);
      ll <- ssm_filter_update_ll(v, Finv);
      K <- ssm_filter_update_K(P, T, Z, Finv);
      // saving
      res[t, 1] <- ll;
      res[t, idx[2, 2]:idx[2, 3]] <- v;
      res[t, idx[3, 2]:idx[3, 3]] <- to_vector_colwise(Finv);
      res[t, idx[4, 2]:idx[4, 3]] <- to_vector_colwise(K);
      res[t, idx[5, 2]:idx[5, 3]] <- a;
      res[t, idx[6, 2]:idx[6, 3]] <- to_vector_colwise(P);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a <- ssm_filter_update_a(a, c, T, v, K);
        P <- ssm_filter_update_P(P, Z, T, Q, R, K);
      }
    }
  }
  return res;
}

int ssm_filter_states_size(int m) {
  int sz;
  sz <- m + m * m;
  return sz;
}

vector ssm_filter_states_get_a(vector x, int m) {
  vector[m] a;
  a <- x[1:m];
  return a;
}

matrix ssm_filter_states_get_P(vector x, int m) {
  matrix[m, m] P;
  P <- to_matrix_colwise(x[(m + 1):(m + m * m)], m, m);
  return P;
}

vector[] ssm_filter_states(vector[] filter, matrix Z) {
  vector[ssm_filter_states_size(cols(Z))] res[size(filter)];
  int n;
  int m;
  int p;
  n <- size(filter);
  m <- cols(Z);
  p <- rows(Z);
  {
    vector[m] aa;
    matrix[m, m] PP;
    vector[p] v;
    matrix[p, p] Finv;
    vector[m] a;
    matrix[m, m] P;

    for (t in 1:n) {
      // updating
      v <- ssm_filter_get_v(filter[t], m, p);
      Finv <- ssm_filter_get_Finv(filter[t], m, p);
      a <- ssm_filter_get_a(filter[t], m, p);
      P <- ssm_filter_get_P(filter[t], m, p);
      aa <- a + P * Z ' * Finv * v;
      PP <- to_symmetric_matrix(P - P * quad_form(Finv, Z) * P);
      // saving
      res[t, 1:m] <- aa;
      res[t, (m + 1):(m + m * m)] <- to_vector_colwise(PP);
    }
  }
  return res;
}


////// Smoothers //////////////////

// ssm_smoother_disturbances
// ssm_smoother_states
// ssm_smoother_states_fast
// ssm_smoother_sim


vector ssm_smooth_update_r(vector r, matrix Z, vector v, matrix Finv,
                           matrix L) {
  vector[num_elements(r)] r_new;
  r_new <- Z' * Finv * v + L ' * r;
  return r_new;
}

matrix ssm_smooth_update_N(matrix N, matrix Z, matrix Finv, matrix L) {
  matrix[rows(N), cols(N)] N_new;
  N_new <- quad_form(Finv, Z) * quad_form(N, L);
  return N_new;
}

int ssm_smooth_state_size(int m) {
  int sz;
  sz <- m + m * m;
  return sz;
}

vector ssm_smooth_state_get_mean(vector x, int m) {
  vector[m] alpha;
  alpha <- x[1:m];
  return alpha;
}

matrix ssm_smooth_state_get_var(vector x, int m) {
  matrix[m, m] V;
  V <- to_matrix_colwise(x[(m + 1):(m + m * m)], m, m);
  return V;
}


// Durbin Koopmans 4.44
vector[] ssm_smooth_state(vector[] filter, matrix Z, matrix T) {
  vector[ssm_smooth_state_size(cols(Z))] res[size(filter)];
  int n;
  int m;
  int p;
  n <- size(filter);
  m <- cols(Z);
  p <- rows(Z);
  {
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[m] alpha;
    matrix[m, m] V;
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    vector[m] a;
    matrix[m, m] P;

    r <- rep_vector(0.0, m);
    N <- rep_matrix(0.0, m, m);
    for (i in 1:n) {
      int t;
      // move backwards in time
      t <- n - i + 1;
      // updating
      K <- ssm_filter_get_K(filter[t], m, p);
      v <- ssm_filter_get_v(filter[t], m, p);
      Finv <- ssm_filter_get_Finv(filter[t], m, p);
      a <- ssm_filter_get_a(filter[t], m, p);
      P <- ssm_filter_get_P(filter[t], m, p);
      L <- ssm_filter_update_L(Z, T, K);
      r <- ssm_smooth_update_r(r, Z, v, Finv, L);
      N <- ssm_smooth_update_N(N, Z, Finv, L);
      alpha <- a + P * r;
      V <- to_symmetric_matrix(P - P * N * P);
      // saving
      res[t, 1:m] <- alpha;
      res[t, (m + 1):(m + m * m)] <- to_vector_colwise(V);
    }
  }
  return res;
}

int ssm_smooth_eps_size(int p) {
  int sz;
  sz <- p + p * p;
  return sz;
}

vector ssm_smooth_eps_get_mean(vector x, int p) {
  vector[p] eps;
  eps <- x[1:p];
  return eps;
}

matrix ssm_smooth_eps_get_var(vector x, int p) {
  matrix[p, p] eps_var;
  eps_var <- to_matrix_colwise(x[(p + 1):(p + p * p)], p, p);
  return eps_var;
}

// Observation disturbance smoother
// Durbin Koopmans Sec 4.5.3 (eq 4.69)
vector[] ssm_smooth_eps(vector[] filter, matrix Z, matrix H, matrix T) {
  vector[ssm_smooth_eps_size(rows(Z))] res[size(filter)];
  int n;
  int m;
  int p;
  n <- size(filter);
  m <- cols(Z);
  p <- rows(Z);
  {
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    vector[p] eps;
    matrix[p, p] var_eps;

    r <- rep_vector(0.0, m);
    N <- rep_matrix(0.0, m, m);
    for (i in 1:n) {
      int t;
      // move backwards in time
      t <- n - i + 1;
      // updating
      K <- ssm_filter_get_K(filter[t], m, p);
      v <- ssm_filter_get_v(filter[t], m, p);
      Finv <- ssm_filter_get_Finv(filter[t], m, p);
      L <- ssm_filter_update_L(Z, T, K);
      r <- ssm_smooth_update_r(r, Z, v, Finv, L);
      N <- ssm_smooth_update_N(N, Z, Finv, L);
      eps <- H * (Finv * v - K ' * r);
      var_eps <- to_symmetric_matrix(H - H * (Finv + quad_form(N, K)) * H);
      // saving
      res[t, 1:p] <- eps;
      res[t, (p + 1):(p + p * p)] <- to_vector_colwise(var_eps);
    }
  }
  return res;
}


int ssm_smooth_eta_size(int q) {
  int sz;
  sz <- q + q * q;
  return sz;
}

vector ssm_smooth_eta_get_mean(vector x, int q) {
  vector[q] eta;
  eta <- x[1:q];
  return eta;
}

matrix ssm_smooth_eta_get_var(vector x, int q) {
  matrix[q, q] eta_var;
  eta_var <- to_matrix_colwise(x[(q + 1):(q + q * q)], q, q);
  return eta_var;
}

// State disturbance smoother
// Durbin Koopmans Sec 4.5.3 (eq 4.69)
vector[] ssm_smooth_eta(vector[] filter,
                        matrix Z, matrix T,
                        matrix R, matrix Q) {
  vector[ssm_smooth_eta_size(cols(Q))] res[size(filter)];
  int n;
  int m;
  int p;
  int q;
  n <- size(filter);
  m <- cols(Z);
  p <- rows(Z);
  q <- rows(Q);
  {
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    vector[q] eta;
    matrix[q, q] var_eta;

    r <- rep_vector(0.0, m);
    N <- rep_matrix(0.0, m, m);
    for (i in 1:n) {
      int t;
      // move backwards in time
      t <- n - i + 1;
      // updating
      K <- ssm_filter_get_K(filter[t], m, p);
      v <- ssm_filter_get_v(filter[t], m, p);
      Finv <- ssm_filter_get_Finv(filter[t], m, p);
      L <- ssm_filter_update_L(Z, T, K);
      r <- ssm_smooth_update_r(r, Z, v, Finv, L);
      N <- ssm_smooth_update_N(N, Z, Finv, L);
      eta <- Q * R ' * r;
      var_eta <- to_symmetric_matrix(Q - Q * quad_form(N, R) * Q);
      // saving
      res[t, 1:q] <- eta;
      res[t, (q + 1):(q + q * q)] <- to_vector_colwise(var_eta);
    }
  }
  return res;
}


// Fast state smoother
// Durbin Koopmans Sec 4.5.3 (eq 4.69)
vector[] ssm_smooth_faststate(vector[] filter,
                              vector c, matrix Z, matrix T,
                              matrix R, matrix Q) {
  vector[cols(Z)] alpha[size(filter)];
  int n;
  int m;
  int p;
  int q;
  n <- size(filter);
  m <- cols(Z);
  p <- rows(Z);
  q <- rows(Q);
  {
    vector[m] r;
    matrix[m, m] L;
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    vector[m] a1;
    matrix[m, m] P1;
    vector[q] eta[n];

    // find smoothed state disturbances
    r <- rep_vector(0.0, m);
    for (i in 1:n) {
      int t;
      // move backwards in time
      t <- n - i + 1;
      // updating
      K <- ssm_filter_get_K(filter[t], m, p);
      v <- ssm_filter_get_v(filter[t], m, p);
      Finv <- ssm_filter_get_Finv(filter[t], m, p);
      L <- ssm_filter_update_L(Z, T, K);
      r <- ssm_smooth_update_r(r, Z, v, Finv, L);
      eta[t] <- Q * R ' * r;
      // saving
    }
    // calculate smoothed states
    a1 <- ssm_filter_get_a(filter[1], m, p);
    P1 <- ssm_filter_get_P(filter[1], m, p);
    alpha[1] <- a1 + P1 * r;
    for (t in 1:(n - 1)) {
      alpha[t + 1] <- c + T * alpha[t] + R * eta[t];
    }
  }
  return alpha;
}

////// Simulators /////////////////
int[,] ssm_sim_idx(int m, int p, int q) {
  int sz[4, 3];
  // y
  sz[1, 1] <- p;
  // a
  sz[2, 1] <- m;
  // eps
  sz[3, 1] <- p;
  // eta
  sz[4, 1] <- q;
  // Fill in start and stop points
  sz[1, 2] <- 1;
  sz[1, 3] <- sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:4) {
    sz[i, 2] <- sz[i - 1, 3] + 1;
    sz[i, 3] <- sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}

int ssm_sim_size(int m, int p, int q) {
  int sz;
  sz <- ssm_sim_idx(m, p, q)[4, 3];
  return sz;
}

vector ssm_sim_get_y(vector x, int m, int p, int q) {
  vector[m] y;
  int idx[4, 3];
  idx <- ssm_sim_idx(m, p, q);
  y <- x[idx[1, 2]:idx[1, 3]];
  return y;
}

vector ssm_sim_get_a(vector x, int m, int p, int q) {
  vector[m] a;
  int idx[4, 3];
  idx <- ssm_sim_idx(m, p, q);
  a <- x[idx[2, 2]:idx[2, 3]];
  return a;
}


vector ssm_sim_get_eps(vector x, int m, int p, int q) {
  vector[m] eps;
  int idx[4, 3];
  idx <- ssm_sim_idx(m, p, q);
  eps <- x[idx[3, 2]:idx[3, 3]];
  return eps;
}

vector ssm_sim_get_eta(vector x, int m, int p, int q) {
  vector[m] eta;
  int idx[4, 3];
  idx <- ssm_sim_idx(m, p, q);
  eta <- x[idx[4, 2]:idx[4, 3]];
  return eta;
}

// only save y and a
vector[] ssm_sim_rng(int n,
                    vector c, matrix Z, matrix H,
                    vector d, matrix T, matrix R, matrix Q,
                    vector a1, matrix P1) {
  vector[ssm_sim_size(cols(Z), rows(Z), cols(Q))] ret[n];
  int p;
  int m;
  int q;
  p <- rows(Z);
  m <- cols(Z);
  q <- cols(Q);
  {
    vector[p] y;
    vector[p] eps;
    vector[m] a;
    vector[q] eta;
    vector[p] zero_p;
    vector[q] zero_q;
    vector[m] zero_m;
    int idx[4, 3];

    idx <- ssm_sim_idx(m, p, q);
    zero_p <- rep_vector(0.0, p);
    zero_q <- rep_vector(0.0, q);
    zero_m <- rep_vector(0.0, m);
    a <- multi_normal_rng(a1, P1);
    for (t in 1:n) {
      eps <- multi_normal_rng(zero_p, H);
      y <- d + Z * a + eps;
      // since eta_t is for alpha_{t + 1}
      if (t == n) {
        eta <- zero_q;
      } else {
        eta <- multi_normal_rng(zero_q, Q);
      }
      // save
      ret[t, idx[1, 2]:idx[1, 3]] <- y;
      ret[t, idx[2, 2]:idx[2, 3]] <- a;
      ret[t, idx[3, 2]:idx[3, 3]] <- eps;
      ret[t, idx[4, 2]:idx[4, 3]] <- eta;
      // a_{t + 1}
      if (t < n) {
        a <- c + T * a + R * eta;
      }
    }
  }
  return ret;
}


// Smoothing Simulators
int[,] ssm_simsmo_dist_idx(int p, int q) {
  int sz[2, 3];
  // eps
  sz[1, 1] <- p;
  // eta
  sz[2, 1] <- q;

  // Fill in start and stop points
  sz[1, 2] <- 1;
  sz[1, 3] <- sz[1, 2] + sz[1, 1] - 1;
  sz[2, 2] <- sz[1, 3] + 1;
  sz[2, 3] <- sz[2, 2] + sz[2, 1] - 1;
  return sz;
}

int ssm_simsmo_dist_size(int p, int q) {
  int sz;
  sz <- ssm_simsmo_dist_idx(p, q)[2, 3];
  return sz;
}

vector ssm_simsmo_get_eta(vector x, int p, int q) {
  int idx[2, 3];
  vector[q] eta;
  idx <- ssm_simsmo_dist_idx(p, q);
  eta <- x[idx[2, 2]:idx[2, 3]];
  return eta;
}

vector ssm_simsmo_get_eps(vector x, int p, int q) {
  int idx[2, 3];
  vector[p] eps;
  idx <- ssm_simsmo_dist_idx(p, q);
  eps <- x[idx[1, 2]:idx[1, 3]];
  return eps;
}

// ssm_simsmo_alpha_rng
vector[] ssm_simsmo_dist_rng(vector[] eps, vector[] eta,
                      vector c, matrix Z, matrix H,
                      vector d, matrix T, matrix R, matrix Q,
                      vector a1, matrix P1) {
    vector[ssm_simsmo_dist_size(rows(T), rows(Q))] draws[size(eps)];
    int n;
    int p;
    int m;
    int q;
    n <- size(eps);
    p <- rows(Z);
    m <- rows(T);
    q <- rows(Q);
    {
      vector[ssm_filter_size(m, p)] filter[n];
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eta_size(p)] epshat_plus[n];
      vector[ssm_smooth_eta_size(q)] etahat_plus[n];
      int idx[2, 3];
      idx <- ssm_simsmo_dist_idx(p, q);
      // simulate unconditional disturbances and observations
      sims <- ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] <- ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter <- ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      epshat_plus <- ssm_smooth_eps(filter, Z, H, T);
      for (i in 1:n) {
        draws[i, idx[1, 2]:idx[1, 3]] <- (ssm_sim_get_eps(sims[i], m, p, q)
                                      - ssm_smooth_eps_get_mean(epshat_plus[i], p)
                                      + ssm_smooth_eps_get_mean(eps[i], p));
      }
      // mean correct eta samples
      etahat_plus <- ssm_smooth_eta(filter, Z, T, R, Q);
      for (i in 1:n) {
        draws[i, idx[2, 2]:idx[2, 3]] <- (ssm_sim_get_eta(sims[i], m, p, q)
                                      - ssm_smooth_eta_get_mean(etahat_plus[i], q)
                                      + ssm_smooth_eta_get_mean(eta[i], q));
      }
    }
    return draws;
}

// ssm_simsmo_alpha_rng
vector[] ssm_simsmo_states_rng(vector[] alpha,
                      vector c, matrix Z, matrix H,
                      vector d, matrix T, matrix R, matrix Q,
                      vector a1, matrix P1) {
    vector[rows(T)] draws[size(alpha)];
    int n;
    int p;
    int m;
    int q;
    n <- size(alpha);
    p <- rows(Z);
    m <- rows(T);
    q <- rows(Q);
    {
      vector[ssm_filter_size(m, p)] filter[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[p] y[n];
      vector[m] alpha_hat_plus[n];
      // simulate unconditional disturbances and observations
      sims <- ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] <- ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter with simulated y's
      filter <- ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      alpha_hat_plus <- ssm_smooth_faststate(filter, c, Z, T, R, Q);
      for (i in 1:n) {
        draws[i] <- (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha[i]);
      }
    }
    return draws;
}

// and de jong methods
// ssm_simsmodj_eta_rng
// ssm_simsmodj_eps_rng
// ssm_simsmodj_alpha_rng

// Partial Autocorrelations to Autocorrelations transformation
// Maps (-1, 1)^p to AR space
// Translated from R stats C function partrans
// https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/pacf.c
// vector pacf_to_acf(vector x) {
//   int n;
//
//   // elements must be between -1 and 1
//   for (i in 1:p) {
//     if (x[i] > 1) {
//       reject("x is greater than 1");
//     } else if (x[i] < -1) {
//       reject("x is less than -1");
//     }
//   }
//
//   work = x;
//   ret = work;
//   /* run the Durbin-Levinson recursions to find phi_{j.},
//      j = 2, ..., p and phi_{p.} are the autoregression coefficients */
//   for(j = 2:p) {
//   	a = ret[j];
//     for(k in 1:j) {
//       // TODO: possible off by one error here
//       work[k] = work[k] - a * ret[j - k - 1];
//     }
//     ret = work;
//   }
//   return ret
// }
//
// // Autocorrelations to Partial Autocorrelations
// // Maps AR space to (-1, 1)^p
// // Translated from R stats C function invpartrans
// // https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/pacf.c
// vector acf_to_pacf(vector x) {
//     real a;
//     vector[size(x)] ret;
//     vector[size(x)] work;
//     int p;
//
//     p = size(x);
//     work = x;
//     for(j in 1:p) {
//       work[j] = phi[j];
//       ret[j] = work[j];
//     }
//     /* Run the Durbin-Levinson recursions backwards
//        to find the PACF phi_{j.} from the autoregression coefficients */
//     for(i in 2:p) {
//       int j;
//       j = p - i;
// 	    a = ret[j];
// 	    for(k in 1:j) {
//         // TODO: possible off by one error here
//         work[k]  = (ret[k] + a * ret[j - k - 1]) / (1 - a * a);
//       }
//       ret = work;
//     }
//     // let user check that outputs are between -1 and 1.
//     return ret;
// }
