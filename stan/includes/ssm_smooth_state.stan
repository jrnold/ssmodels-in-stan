///
/// # State Smoother
///

/**
The number of elements in vectors returned by `ssm_smooth_state`

@param int m The number of states.
@return int The size of the vectors is $m + m ^ 2$.

*/
int ssm_smooth_state_size(int m) {
  int sz;
  sz = m + m * m;
  return sz;
}

/**
Extract $\hat{\vec{\alpha}}_t$ from vectors returned by `ssm_smooth_state`

@param vector x A vector returned by `ssm_smooth_state`
@param int q The number of state disturbances, $\vec{\eta}_t$.
@return vector An $m \times 1$ vector with $\hat{\vec{\eta}}_t$.
*/
vector ssm_smooth_state_get_mean(vector x, int m) {
  vector[m] alpha;
  alpha = x[1:m];
  return alpha;
}

/**
Extract $mat{V}_t$ from vectors returned by `ssm_smooth_state`

@param vector x A vector returned by `ssm_smooth_state`
@param int m The number of states
@return matrix An $m \times m$ matrix with $\mat{V}_t$.
*/
matrix ssm_smooth_state_get_var(vector x, int m) {
  matrix[m, m] V;
  V = to_matrix_colwise(x[(m + 1):(m + m * m)], m, m);
  return V;
}


/**
The state smoother

This calculates the mean and variance of the states, $\vec{\alpha}_t$, given the entire sequence, $\vec{y}_{1:n}$.

@param vector[] filter Results of `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@return vector[] An array of vectors constaining $\hat{\vec{\alpha}}_t$ and $\mat{V}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:n})$.
  in the format described below.

For `Z` and `T` the array can have a size of 1, if it is not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `T`) if it is time varying.

The vectors returned by this function have $m + m ^ 2$ elements in this format,
$$
(\hat{\vec{\alpha}}_t', \VEC\mat{V}_t' ).
$$
Use the `ssm_smooth_state_get_mean` and `ssm_smooth_state_get_var` to extract components
from the returned vectors.

See [@DurbinKoopman2012, Eq 4.44 (eq 4.69)]
*/
vector[] ssm_smooth_state(vector[] filter, matrix[] Z, matrix[] T) {
  vector[ssm_smooth_state_size(dims(Z)[3])] res[size(filter)];
  int n;
  int m;
  int p;
  n = size(filter);
  m = dims(Z)[3];
  p = dims(Z)[2];
  {
    // system matrices for current iteration
    matrix[p, m] Z_t;
    matrix[m, m] T_t;
    // smoother results
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[m] alpha;
    matrix[m, m] V;
    // filter results
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    vector[m] a;
    matrix[m, m] P;

    if (size(Z) == 1) {
      Z_t = Z[1];
    }
    if (size(T) == 1) {
      T_t = T[1];
    }
    // initialize smoother
    r = rep_vector(0.0, m);
    N = rep_matrix(0.0, m, m);
    for (i in 0:(n - 1)) {
      int t;
      // move backwards in time
      t = n - i;
      // set time-varying system matrices
      if (size(Z) > 1) {
        Z_t = Z[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      // get filtered values
      K = ssm_filter_get_K(filter[t], m, p);
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      a = ssm_filter_get_a(filter[t], m, p);
      P = ssm_filter_get_P(filter[t], m, p);
      // updating
      L = ssm_filter_update_L(Z_t, T_t, K);
      r = ssm_smooth_update_r(r, Z_t, v, Finv, L);
      N = ssm_smooth_update_N(N, Z_t, Finv, L);
      alpha = a + P * r;
      V = to_symmetric_matrix(P - P * N * P);
      // saving
      res[t, 1:m] = alpha;
      res[t, (m + 1):(m + m * m)] = to_vector(V);
    }
  }
  return res;
}
