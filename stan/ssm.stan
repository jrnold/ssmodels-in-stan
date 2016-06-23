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
int[,] ssm_filter_return_idx(int m, int p) {
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

int ssm_filter_return_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx <- ssm_filter_return_idx(m, p);
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

real ssm_loglik(vector[] y,
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


// Filtering
vector[] ssm_filter(
                vector[] y,
                vector c, matrix Z, matrix H,
                vector d, matrix T, matrix R, matrix Q,
                vector a1, matrix P1) {

  // returned data
  vector[ssm_filter_return_size(cols(Z), rows(Z))] res[size(y)];
  int r;
  int n;
  int p;
  int m;

  // sizes
  n <- size(y); // number of obs
  p <- rows(Z); // obs size
  m <- cols(Z); // number of states
  r <- cols(Q); // number of state disturbances

  //print("Sizes: n = ", m, ", p = ", n, ", m = ", m, ", r = ", r);
  {
    vector[m] a;
    matrix[m, m] P;
    vector[m * m] P_vec;
    vector[p] v;
    matrix[p, p] Finv;
    vector[p * p] Finv_vec;
    matrix[m, p] K;
    vector[m * p] K_vec;
    real ll;
    int idx[6, 3];

    idx <- ssm_filter_return_idx(m, p);
    a <- a1;
    res[1, idx[5, 2]:idx[5, 3]] <- a;
    P <- P1;
    P_vec <- to_vector_colwise(P);
    res[1, idx[6, 2]:idx[6, 3]] <- P_vec;
    for (t in 1:n) {
      v <- ssm_filter_update_v(y[t], a, d, Z);
      res[t, idx[2, 2]:idx[2, 3]] <- v;
      Finv <- ssm_filter_update_Finv(P, Z, H);
      Finv_vec <- to_vector_colwise(Finv);
      res[t, idx[3, 2]:idx[3, 3]] <- Finv_vec;
      K <- ssm_filter_update_K(P, T, Z, Finv);
      K_vec <- to_vector_colwise(K);
      res[t, idx[4, 2]:idx[4, 3]] <- K_vec;
      ll <- ssm_filter_update_ll(v, Finv);
      res[t, 1] <- ll;
      // don't save a, P for last iteration
      if (t < n) {
        a <- ssm_filter_update_a(a, c, T, v, K);
        res[t + 1, idx[5, 2]:idx[5, 3]] <- a;
        P <- ssm_filter_update_P(P, Z, T, Q, R, K);
        P_vec <- to_vector_colwise(P);
        res[t + 1, idx[6, 2]:idx[6, 3]] <- P_vec;
      }
    }
  }
  return res;
}

////// Smoothers //////////////////

// ssm_smoother_disturbances
// ssm_smoother_states
// ssm_smoother_states_fast
// ssm_smoother_sim

////// Simulators /////////////////
int[,] ssm_sim_return_idx(int m, int p) {
  int sz[6, 3];
  // y
  sz[1, 1] <- p;
  // a
  sz[2, 1] <- m;
  // Fill in start and stop points
  sz[1, 2] <- 1;
  sz[1, 3] <- sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:2) {
    sz[i, 2] <- sz[i - 1, 3] + 1;
    sz[i, 3] <- sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}

int ssm_sim_return_size(int m, int p) {
  int sz;
  sz <- ssm_sim_return_idx(m, p)[2, 3];
  return sz;
}

// only save y and a
vector[] ssm_sim_rng(int n,
                    vector c, matrix Z, matrix H,
                    vector d, matrix T, matrix R, matrix Q,
                    vector a1, matrix P1) {
  vector[ssm_sim_return_size(cols(Z), rows(Z))] ret[n];
  int p;
  int m;
  int r;
  p <- rows(Z);
  m <- cols(Z);
  r <- cols(Q);
  {
    vector[p] y;
    vector[p] eps;
    vector[m] a;
    vector[m] eta;
    matrix[m, m] RQR;
    int idx[2, 3];

    idx <- ssm_sim_return_idx(m, p);
    RQR <- quad_form_sym(Q, R);
    a <- multi_normal_rng(a1, P1);
    ret[1, idx[2, 2]:idx[2, 3]] <- a;
    for (t in 1:n) {
      y <- multi_normal_rng(d + Z * a, H);
      ret[1, idx[1, 2]:idx[1, 3]] <- y;
      if (t < n) {
        a <- multi_normal_rng(c + T * a, RQR);
        ret[1, idx[2, 2]:idx[2, 3]] <- a;
      }
    }
  }
  return ret;
}

// simulate from fixed starting values
vector[] ssm_sim2_rng(int n,
                    vector c, matrix Z, matrix H,
                    vector d, matrix T, matrix R, matrix Q,
                    vector a1) {
  vector[ssm_sim_return_size(cols(Z), rows(Z))] ret[n];
  int p;
  int m;
  int r;
  p <- rows(Z);
  m <- cols(Z);
  r <- cols(Q);
  {
    vector[p] y;
    vector[p] eps;
    vector[m] a;
    vector[m] eta;
    matrix[m, m] RQR;
    int idx[2, 3];

    idx <- ssm_sim_return_idx(m, p);
    RQR <- quad_form_sym(Q, R);
    a <- a1;
    ret[1, idx[2, 2]:idx[2, 3]] <- a;
    for (t in 1:n) {
      y <- multi_normal_rng(d + Z * a, H);
      ret[1, idx[1, 2]:idx[1, 3]] <- y;
      if (t < n) {
        a <- multi_normal_rng(c + T * a, RQR);
        ret[1, idx[2, 2]:idx[2, 3]] <- a;
      }
    }
  }
  return ret;
}


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
