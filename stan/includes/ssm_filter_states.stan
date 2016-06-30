///
/// # Get Filtered states
///
/// Requires: `ssm_utils.stan`, `ssm_filter.stan`

/**
Length of the vectors returned by `ssm_filter_states`

@param int m Number of states
@return int The size of the vector

*/
int ssm_filter_states_size(int m) {
  int sz;
  sz = m + m * m;
  return sz;
}

/**
Extract $a_{t|t}$ from the results of `ssm_filter_states`

@param vector x A vector returned by `ssm_filter_states`
@param int m Number of states
@return matrix An $m \times 1$ vector with the filtered expected value of the staate,
  $\vec{a}_{t|t} = \E(\vec{\alpha}_t | \vec{y}_{1:t})$.

*/
vector ssm_filter_states_get_a(vector x, int m) {
  vector[m] a;
  a = x[1:m];
  return a;
}

/**
Extract $P_{t|t}$ from the results of `ssm_filter_states`

@param vector x A vector returned by `ssm_filter_states`
@param int m Number of states
@return matrix An $m \times m$ matrix with the filtered variance of the staate,
  $\mat{P}_{t|t} = \Var(\vec{\alpha}_t | \vec{y}_{1:t})$.

*/
matrix ssm_filter_states_get_P(vector x, int m) {
  matrix[m, m] P;
  P = to_matrix_colwise(x[(m + 1):(m + m * m)], m, m);
  return P;
}

/**
Calculate filtered expected values and variances of the states

The filtering function `ssm_filter` returns the mean and variance of the predicted states,
$\vec{a}_t = \E(\vec{alpha}_t | \vec{y}_{1:(t -1)})$ and $\mat{P}_t = \Var(\vec{alpha}_t | \vec{y}_{1:(t -1)})$.
Using these results, this function calculates the mean and variance of the filtered states, $\vec{a}_{t|t} = \E(\vec{alpha}_t | \vec{y}_{1:t})$ and $\mat{P}_{t|t} = \Var(\vec{alpha}_t | \vec{y}_{1:t})$.

@param vector[] filter Results from `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@return Array of size $n$ of vectors.

The vectors returned by `ssm_filter_states` are of length $m + m ^ 2$, with
$$
v_t = (\vec{a}_{t|t}', \VEC\vec{P}_{t|t})'
$$
Use the functions `ssm_filter_states_get_a` and `ssm_filter_states_get_P` to extract
elements from the results.

For `Z` the array can have a size of 1, if it is not time-varying, or a size of $n - 1$ if it is time varying.

*/
vector[] ssm_filter_states(vector[] filter, matrix[] Z) {
  vector[ssm_filter_states_size(dims(Z)[3])] res[size(filter)];
  int n;
  int m;
  int p;
  n = size(filter);
  m = dims(Z)[3];
  p = dims(Z)[2];
  {
    // system matrices for current iteration
    matrix[p, m] Z_t;
    // filter matrices
    vector[m] aa; // filtered values of the state, a_{t|t}
    matrix[m, m] PP; // filtered values of the variance of the state, P_{t|t}
    vector[p] v;
    matrix[p, p] Finv;
    vector[m] a;
    matrix[m, m] P;

    Z_t = Z[1];
    for (t in 1:n) {
      if (t > 1) {
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
      }
      // extract values from the filter
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      a = ssm_filter_get_a(filter[t], m, p);
      P = ssm_filter_get_P(filter[t], m, p);
      // calcualte filtered values
      aa = a + P * Z_t ' * Finv * v;
      PP = to_symmetric_matrix(P - P * quad_form(Finv, Z_t) * P);
      // saving
      res[t, 1:m] = aa;
      res[t, (m + 1):(m + m * m)] = to_vector(PP);
    }
  }
  return res;
}
