///
/// # Filtering
///
/// **Requires**: `ssm_utils.stan`, `ssm_filter_utils.stan`.

/**
Indexes of the return values of the Kalman filter functions:
`ssm_filter`.

@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return int[,] A $6 \times 3$ integer array containing the indexes of the return values of the Kalman filter.

`ssm_filter_idx` returns a $6 \times 3$ integer array with the
(length, start index, stop index) of ($L$, $\vec{v}$, $\vec{F}^-1$, $\mat{K}$, $\vec{a}$, $\mat{P}$).

value    length    start                 stop
-------- --------- --------------------- ----------------------
log-lik  1         1                     1
v        p         2                     1 + p
F^-1     p^2       2 + p                 1 + p + p ^ 2
K        mp        2 + p + p^2           1 + p + p ^ 2 + mp
a_t      m         2 + p + p^2 + mp      1 + p + p ^ 2 + mp + m
P^t      m * m     2 + p + p^2 + mp + m  1 + p + p ^ 2 + mp + m ^ 2

*/
int[,] ssm_filter_idx(int m, int p) {
  int sz[6, 3];
  // loglike
  sz[1, 1] = 1;
  // v
  sz[2, 1] = p;
  // Finv
  sz[3, 1] = p * p;
  // K
  sz[4, 1] = m * p;
  // a
  sz[5, 1] = m;
  // P
  sz[6, 1] = m * m;
  // Fill in start and stop points
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:6) {
    sz[i, 2] = sz[i - 1, 3] + 1;
    sz[i, 3] = sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}

/**
Number of elements in vector containing filter results

@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return int The number of elements in the vector.

*/
int ssm_filter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  sz = idx[6, 3];
  return sz;
}

/**
Get the log-likehood from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return real The log-likelihood $L_t$

*/
real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y = x[1];
  return y;
}

/**
Get the forecast error from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return vector A $p \times 1$ vector with the forecast error, $v_t$.

*/
vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[2, 2], idx[2, 3]);
  return y;
}

/**
Get the forecast precision from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return matrix A $p \times p$ matrix with the forecast precision, $\mat{F}^{-1}_t$.

*/
matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[3, 2], idx[3, 3]), p, p);
  return y;
}

/**
Get the Kalman gain from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return matrix A $m \times p$ matrix with the Kalman gain, $\mat{F}^{-1}_t$.

*/
matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[4, 2], idx[4, 3]), m, p);
  return y;
}

/**
Get the expected value of the predicted state from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return vector An $m \times 1$ vector with the expected value of the predicted state, $\E(\vec{alpha}_t | \vec{y}_{1:(t-1)} = a_t$.

*/
vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[5, 2], idx[5, 3]);
  return y;
}

/**
Get the variance of the predicted state from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return matrix An $m \times m$ matrix with the variance of the predicted state, $\Var(\vec{alpha}_t | \vec{y}_{1:(t-1)} = P_t$.

*/
matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[6, 2], idx[6, 3]), m, m);
  return y;
}

/**
Kalman filter

@param vector[] y Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
@param vector[] d Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param vector[] c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@param vector a1 Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
@param matrix P1 Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
@return vector[] Array of size $n$ of $(1 + p + p^2 + mp + m + m^2) \times 1$ vectors in the format described in `ssm_filter_idx`.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

`ssm_filter` runs a forward filter on the state space model and calculates,

- log-likelihood ($\ell_t$) for each observation
- Forecast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y}_{1:(t -1)})$.
- Forecast precision, $\mat{F}^{-1}_t$.
- Kalman gain, $\mat{K}_t$.
- Predicted states, $a_t = \E(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.
- Variance of the predicted states, $P_t = \Var(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.

The results of Kalman filter for a given are returned as a $1 + p + p ^ 2 + mp + m ^ 2$ vector for each time period, where
$$
(\ell_t, \vec{v}_t', \VEC\mat{F}^{-1}_t', \VEC\mat{K}_t', \vec{a}_t', \VEC\mat{P}_t' )''.
$$

*/
vector[] ssm_filter(vector[] y,
                    vector[] d, matrix[] Z, matrix[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {

  // returned data
  vector[ssm_filter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
  int q;
  int n;
  int p;
  int m;

  // sizes
  n = size(y); // number of obs
  p = dims(Z)[2]; // obs size
  m = dims(Z)[3]; // number of states
  q = dims(Q)[2]; // number of state disturbances

  //print("Sizes: n = ", m, ", p = ", n, ", m = ", m, ", q = ", q);
  {
    // system matrices for current iteration
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    // result matricees for each iteration
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    real ll;
    int idx[6, 3];

    idx = ssm_filter_idx(m, p);

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form(Q_t, R_t);
    a = a1;
    P = P1;
    for (t in 1:n) {
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
        }
        if (size(c) > 1) {
          c_t = c[t];
        }
        if (size(T) > 1) {
          T_t = T[t];
        }
        if (size(R) > 1) {
          R_t = R[t];
        }
        if (size(Q) > 1) {
          Q_t = Q[t];
        }
        if (size(R) > 1 && size(Q) > 1) {
          RQR = quad_form(Q_t, R_t);
        }
      }
      // updating
      v = ssm_filter_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_filter_update_Finv(P, Z_t, H_t);
      K = ssm_filter_update_K(P, T_t, Z_t, Finv);
      ll = ssm_filter_update_ll(v, Finv);
      // saving
      res[t, 1] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = to_vector(Finv);
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = to_vector(P);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a = ssm_filter_update_a(a, c_t, T_t, v, K);
        P = ssm_filter_update_P(P, Z_t, T_t, RQR, K);
      }
    }
  }
  return res;
}
