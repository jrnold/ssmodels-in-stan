// -*- mode: stan -*-

// Utility Functions

matrix make_symmetric_matrix(matrix x) {
  return 0.5 * (x + x');
}

matrix to_matrix_colwise(vector v, int r, int c) {
  matrix[r, c] res;
  for (i in 1:size(v)) {
      res[(i - 1) % r + 1, (i - 1) / r + 1] <- v[i]
  }
  return res;
}

matrix to_matrix_rowwise(vector v, int r, int c) {
  matrix[r, c] res;
  for (i in 0:(size(v) - 1)) {
      res[(i - 1) / c + 1, (i - 1) % c + 1] <- v[i]
  }
  return res;
}

matrix to_vector_colwise(matrix x) {
  vector[num_elements(x)] res;
  for (r in 1:rows(x)) {
    for (c in 1:cols(x)) {
      res[r * (c - 1) + r] <- x[r, c]
    }
  }
  return res;
}

matrix to_vector_rowwise(matrix x) {
  vector[num_elements(x)] res;
  for (r in 1:rows(x)) {
    for (c in 1:cols(x)) {
      res[(r - 1) * c + c] <- x[r, c]
    }
  }
  return res;
}

// Length of vectors that SSM returns
// value    size      location
// log-lik  1         1
// v        p         2
// F^-1     p^2       2 + p
// K        mp        2 + p + p^2
// a_t      m         2 + p + p^2 + mp
// P^t      m         2 + p + p^2 + mp + m
int ssm_filter_return_size(int m, int p) {
  int sz;
  sz <- 1 + p + p * p + m + m * m + m * p
  return sz;
}

real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y <- x[1];
  return y;
}

vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] v;
  y <- segment(x, 2, p);
  return y;
}

matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p * p] y;
  y <- to_matrix_colwise(segment(x, 2 + p, p * p), p, p);
  return y;
}

matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[p * p] y;
  y <- to_matrix_colwise(segment(x, 2 + p + p * p, m * p), m, p);
  return y;
}

vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  y <- segment(x, 2 + p + p * p + m * p, m);
  return y;
}

matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m * m] y;
  y <- to_matrix_colwise(segment(x,  2 + p + p * p + m * p + m, m * m), m, m);
  return y;
}

// Filtering
vector[] ssm_filter(vector[] y,
                    vector c, matrix Z, matrix H,
                    vector d, matrix T, matrix R, matrix Q,
                    vector a1, matrix P1) {

  // p = rows(Z)
  // r = cols(Z)
  // n = size(y)
  // returned data
  vector[ssf_filter_return_size(rows(Z), cols(Z))] res[size(y)];
  // temp variables
  int p;
  int r;
  int n;
  vector[rows(Z)] v;
  matrix[rows(Z), rows(Z)] F;
  matrix[rows(Z), rows(Z)] Finv;
  matrix[rows(Z), cols(Z)] K;

  // sizes
  n <- size(y); // number of obs
  p <- rows(Z); // obs size
  m <- cols(Z); // number of states

  for (t in 1:n) {
    // PREDICT STATES
    // one step ahead predictive distribion of p(\theta_t | y_{1:(t-1)})
    a <- g[t] + G[t] * m;
    R <- quad_form(C, G[t] ') + W[t];
    m <- a;
    C <- R;
    // print("t=", t, ", a=", a);
    // print("t=", t, ", R=", R);
    for (j in 1:r) {
      if (int_step(miss[t, j])) {
        e[j] <- 0.0;
        Q_inv[j] <- 0.0;
      } else {
        // print("t = ", t, ", predict");
        // PREDICT OBS
        // one step ahead predictive distribion of p(y_t | y_{1:(t-1)})
        Fj <- row(F[t], j) ';
        f <- b[t, j] + dot_product(Fj, m);
        Q <- quad_form(C, Fj) + V[t, j];
        Q_inv[j] <- 1.0 / Q;
        // forecast error
        e[j] <- y[t, j] - f;
        // Kalman gain
        K <- C * Fj * Q_inv[j];
        // print("t = ", t, ", filter");
        // FILTER STATES
        // posterior distribution of p(\theta_t | y_{1:t})
        m <- m + K * e[j];
        C <- make_symmetric_matrix(C - K * Q * K ');
      }
      // print("t=", t, ", j=", j, ", m=", m);
      // print("t=", t, ", j=", j, ", C=", C);
      // print("t=", t, ", j=", j, ", Q_inv=", Q_inv);
      // print("t=", t, ", j=", j, ", f=", f);
      // print(" ");
    }
    for(i in 1:p) {
      res[t + 1, i] <- m[i];
    }
    C_vec <- to_vector(C);
    for (i in 1:(p * p)) {
      res[t + 1, p + i] <- C_vec[i];
    }
    for(i in 1:p) {
      res[t + 1, p + p * p + i] <- a[i];
    }
    R_vec <- to_vector(R);
    for (i in 1:(p * p)) {
      res[t + 1, 2 * p + p * p + i] <- R_vec[i];
    }
    for (i in 1:r) {
      res[t + 1, 2 * p + 2 * p * p + i] <- e[i];
    }
    for (i in 1:r) {
      res[t + 1, 2 * p + 2 * p * p + r + i] <- Q_inv[i];
    }
  }
  return res;
}

// // Partial Autocorrelations to Autocorrelations transformation
// // Maps (-1, 1)^p to AR space
// // Translated from R stats C function partrans
// // https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/pacf.c
// vector pacf_to_acf(vector x) {
//   vector[size(x)] ret;
//   vector[size(x)] work;
//   real a;
//   int p;
//   p = size(x);
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
