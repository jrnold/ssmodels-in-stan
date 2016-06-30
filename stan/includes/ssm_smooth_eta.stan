/**
The size of the vectors returned by `ssm_smooth_eta`

@param int p The length of the observation vectors, $\vec{y}_t$.
@return int The size of the vectors is $q + q ^ 2$.

*/
int ssm_smooth_eta_size(int q) {
  int sz;
  sz = q + q * q;
  return sz;
}

/**
Extract $\hat{\vec{\varepsilon}}_t$ from vectors returned by `ssm_smooth_eta`

@param vector x A vector returned by `ssm_smooth_eta`
@param int q The number of state disturbances, $\vec{\eta}_t$.
@return vector A $q \times 1$ vector with $\hat{\vec{\eta}}_t$.
*/
vector ssm_smooth_eta_get_mean(vector x, int q) {
  vector[q] eta;
  eta = x[1:q];
  return eta;
}

/**
Extract $\Var(\eta_t|\vec{y}_{1:n})$ from vectors returned by `ssm_smooth_eta`

@param vector x A vector returned by `ssm_smooth_eta`
@param int q The number of state disturbances, $\vec{\eta}_t$.
@return matrix A $q \times q$ matrix with $\Var{\vec{\eta}_t | \vec{y}_{1:n}}$.

*/
matrix ssm_smooth_eta_get_var(vector x, int q) {
  matrix[q, q] eta_var;
  eta_var = to_matrix_colwise(x[(q + 1):(q + q * q)], q, q);
  return eta_var;
}

/**
The state disturbance smoother

This calculates the mean and variance of the observation disturbances, $\vec{\eta}_t$,
given the entire sequence, $\vec{y}_{1:n}$.

@param vector[] filter Results of `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@return vector[] An array of vectors constaining $\hat{\vec{eta}}_t$ and $\Var(\vec{\eta}_t | \vec{y}_{1:n})$
  in the format described below.

For `Z`, `T`, `R`, `Q` the array can have a size of 1, if it is not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `T`, `R`, `Q`) if it is time varying.

The vectors returned by this function have $q + q ^ 2$ elements in this format,
$$
(\hat{\vec{eta}}_t', \VEC\Var(\vec{\eta}_t | \vec{y}_{1:n})' ).
$$
Use the `ssm_smooth_eta_get_mean` and `ssm_smooth_eta_get_var` to extract components
from the returned vectors.

See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]

*/
vector[] ssm_smooth_eta(vector[] filter,
                        matrix[] Z, matrix[] T,
                        matrix[] R, matrix[] Q) {
  vector[ssm_smooth_eta_size(dims(Q)[2])] res[size(filter)];
  int n;
  int m;
  int p;
  int q;
  n = size(filter);
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    // smoother matrices
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[q] eta;
    matrix[q, q] var_eta;
    // system matrices
    matrix[p, m] Z_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    // filter matrices
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;

    // set time-invariant matrices
    if (size(Z) == 1) {
      Z_t = Z[1];
    }
    if (size(T) == 1) {
      T_t = T[1];
    }
    if (size(R) == 1) {
      R_t = R[1];
    }
    if (size(Q) == 1) {
      Q_t = Q[1];
    }
    // initialize smoother
    r = rep_vector(0.0, m);
    N = rep_matrix(0.0, m, m);
    for (i in 1:n) {
      int t;
      // move backwards in time
      t = n - i + 1;
      // update time-varying system matrices
      if (size(Z) > 1) {
        Z_t = Z[t];
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
      // get values from filter
      K = ssm_filter_get_K(filter[t], m, p);
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      // update smoother
      L = ssm_filter_update_L(Z_t, T_t, K);
      r = ssm_smooth_update_r(r, Z_t, v, Finv, L);
      N = ssm_smooth_update_N(N, Z_t, Finv, L);
      eta = Q_t * R_t ' * r;
      var_eta = to_symmetric_matrix(Q_t - Q_t * quad_form(N, R_t) * Q_t);
      // saving
      res[t, 1:q] = eta;
      res[t, (q + 1):(q + q * q)] = to_vector(var_eta);
    }
  }
  return res;
}
