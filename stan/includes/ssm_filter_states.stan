///
/// # Get Filtered states
///
/// Requires: `ssm_utils.stan`, `ssm_filter.stan`

int ssm_filter_states_size(int m) {
  int sz;
  sz = m + m * m;
  return sz;
}

vector ssm_filter_states_get_a(vector x, int m) {
  vector[m] a;
  a = x[1:m];
  return a;
}

matrix ssm_filter_states_get_P(vector x, int m) {
  matrix[m, m] P;
  P = to_matrix_colwise(x[(m + 1):(m + m * m)], m, m);
  return P;
}

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
      // updating
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      a = ssm_filter_get_a(filter[t], m, p);
      P = ssm_filter_get_P(filter[t], m, p);
      aa = a + P * Z_t ' * Finv * v;
      PP = to_symmetric_matrix(P - P * quad_form(Finv, Z_t) * P);
      // saving
      res[t, 1:m] = aa;
      res[t, (m + 1):(m + m * m)] = to_vector(PP);
    }
  }
  return res;
}
