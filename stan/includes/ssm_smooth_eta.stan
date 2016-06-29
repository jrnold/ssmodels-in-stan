int ssm_smooth_eta_size(int q) {
  int sz;
  sz = q + q * q;
  return sz;
}

vector ssm_smooth_eta_get_mean(vector x, int q) {
  vector[q] eta;
  eta = x[1:q];
  return eta;
}

matrix ssm_smooth_eta_get_var(vector x, int q) {
  matrix[q, q] eta_var;
  eta_var = to_matrix_colwise(x[(q + 1):(q + q * q)], q, q);
  return eta_var;
}

/** State disturbance smoother

See [@DurbinKoopman2012, Sec 4.5.3, eq 4.69]
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
