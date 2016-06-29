/** Fast state smoother

The fast state smoother calculates $\hat{\vec{alpha}}_t = \E(\vec{alpha}_t | \vec{y}_{1:n})$.
$$
\hat{\vec{alpha}}_{t + 1} = T_t \hat{\vec{alpha}}_{t} + \mat{R}_t \mat{Q}_t \mat{R}'_t \vec{r}_t ,
$$
where $r_t$ is calcualted from the state disturbance smoother.
The smoother is initialized at $t = 1$ with $\hat{\vec{\alpha}}_t = \vec{a}_1 + \mat{P}_1 \vec{r}_0$.

Unlike the normal state smoother, it does not calculate the variances of the smoothed state.

See [@DurbinKoopmans2012, Sec 4.5.3 (eq 4.69)]

*/
vector[] ssm_smooth_faststate(vector[] filter,
                              vector[] c, matrix[] Z, matrix[] T,
                              matrix[] R, matrix[] Q) {
  vector[dims(Z)[3]] alpha[size(filter)];
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
    vector[m] r[n];
    matrix[m, m] L;
    vector[m] a1;
    matrix[m, m] P1;
    // filter matrices
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    // system matrices
    matrix[p, m] Z_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[p, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    // set time-invariant matrices
    if (size(c) == 1) {
      c_t = c[1];
    }
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
    if (size(Q) == 1 && size(R) == 1) {
      RQR = quad_form(Q[1], R[1]');
    }
    // find smoothed state disturbances
    // Since I don't need to calculate the
    // variances of the smoothed disturbances,
    // I reimplement the state distrurbance smoother here
    // removing extraneous parts.
    r = rep_vector(0.0, m);
    for (i in 0:(n - 1)) {
      int t;
      // move backwards in time
      t = n - i;
      // update time varying system matrices
      if (size(Z) > 1) {
        Z_t = Z[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      // get filter values
      K = ssm_filter_get_K(filter[t], m, p);
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      // updating smoother
      L = ssm_filter_update_L(Z_t, T_t, K);
      r = ssm_smooth_update_r(r, Z_t, v, Finv, L);
    }
    // calculate smoothed states
    a1 = ssm_filter_get_a(filter[1], m, p);
    P1 = ssm_filter_get_P(filter[1], m, p);
    alpha[1] = a1 + P1 * r[1];
    for (t in 1:(n - 1)) {
      if (size(c) > 1) {
        c_t = c[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      if (size(Q) > 1) {
        Q_t = Q[t];
      }
      if (size(R) > 1) {
        R_t = R[t];
      }
      if (size(Q) > 1 || size(R) > 1) {
        RQR = quad_form(Q_t, R_t');
      }
      alpha[t + 1] = c_t + T_t * alpha[t] + RQR * r[t];
    }
  }
  return alpha;
}
