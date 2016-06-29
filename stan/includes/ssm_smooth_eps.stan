int ssm_smooth_eps_size(int p) {
  int sz;
  sz = p + p * p;
  return sz;
}

vector ssm_smooth_eps_get_mean(vector x, int p) {
  vector[p] eps;
  eps = x[1:p];
  return eps;
}

/**

*/
matrix ssm_smooth_eps_get_var(vector x, int p) {
  matrix[p, p] eps_var;
  eps_var = to_matrix_colwise(x[(p + 1):(p + p * p)], p, p);
  return eps_var;
}

/**

Observation disturbance smoother
See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]

*/
vector[] ssm_smooth_eps(vector[] filter, matrix[] Z, matrix[] H, matrix[] T) {
  vector[ssm_smooth_eps_size(dims(Z)[2])] res[size(filter)];
  int n;
  int m;
  int p;
  n = size(filter);
  m = dims(Z)[3];
  p = dims(Z)[2];
  {
    // smoother values
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[p] eps;
    matrix[p, p] var_eps;
    // filter results
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    // system matrices
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    matrix[m, m] T_t;

    // set matrices if time-invariant
    if (size(Z) == 1) {
      Z_t = Z[1];
    }
    if (size(H) == 1) {
      H_t = H[1];
    }
    if (size(T) == 1) {
      T_t = T[1];
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
      if (size(H) > 1) {
        H_t = H[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      // get values from filter
      K = ssm_filter_get_K(filter[t], m, p);
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      // updating
      L = ssm_filter_update_L(Z_t, T_t, K);
      r = ssm_smooth_update_r(r, Z_t, v, Finv, L);
      N = ssm_smooth_update_N(N, Z_t, Finv, L);
      eps = H_t * (Finv * v - K ' * r);
      var_eps = to_symmetric_matrix(H_t - H_t * (Finv + quad_form(N, K)) * H_t);
      // saving
      res[t, 1:p] = eps;
      res[t, (p + 1):(p + p * p)] = to_vector(var_eps);
    }
  }
  return res;
}
