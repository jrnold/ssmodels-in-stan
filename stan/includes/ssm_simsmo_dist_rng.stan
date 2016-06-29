// Smoothing Simulators
int[,] ssm_simsmo_dist_idx(int p, int q) {
  int sz[2, 3];
  // eps
  sz[1, 1] = p;
  // eta
  sz[2, 1] = q;

  // Fill in start and stop points
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  sz[2, 2] = sz[1, 3] + 1;
  sz[2, 3] = sz[2, 2] + sz[2, 1] - 1;
  return sz;
}

int ssm_simsmo_dist_size(int p, int q) {
  int sz;
  sz = ssm_simsmo_dist_idx(p, q)[2, 3];
  return sz;
}

vector ssm_simsmo_get_eta(vector x, int p, int q) {
  int idx[2, 3];
  vector[q] eta;
  idx = ssm_simsmo_dist_idx(p, q);
  eta = x[idx[2, 2]:idx[2, 3]];
  return eta;
}

vector ssm_simsmo_get_eps(vector x, int p, int q) {
  int idx[2, 3];
  vector[p] eps;
  idx = ssm_simsmo_dist_idx(p, q);
  eps = x[idx[1, 2]:idx[1, 3]];
  return eps;
}

// ssm_simsmo_alpha_rng
vector[] ssm_simsmo_dist_rng(vector[] eps, vector[] eta,
                      vector[] d, matrix[] Z, matrix[] H,
                      vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                      vector a1, matrix P1) {
    vector[ssm_simsmo_dist_size(dims(Z)[3], dims(Q)[2])] draws[size(eps)];
    int n;
    int p;
    int m;
    int q;
    n = size(eps);
    p = dims(Z)[2];
    m = dims(Z)[3];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter[n];
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eta_size(p)] epshat_plus[n];
      vector[ssm_smooth_eta_size(q)] etahat_plus[n];
      int idx[2, 3];
      idx = ssm_simsmo_dist_idx(p, q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      epshat_plus = ssm_smooth_eps(filter, Z, H, T);
      for (i in 1:n) {
        draws[i, idx[1, 2]:idx[1, 3]] = (ssm_sim_get_eps(sims[i], m, p, q)
                                      - ssm_smooth_eps_get_mean(epshat_plus[i], p)
                                      + ssm_smooth_eps_get_mean(eps[i], p));
      }
      // mean correct eta samples
      etahat_plus = ssm_smooth_eta(filter, Z, T, R, Q);
      for (i in 1:n) {
        draws[i, idx[2, 2]:idx[2, 3]] = (ssm_sim_get_eta(sims[i], m, p, q)
                                      - ssm_smooth_eta_get_mean(etahat_plus[i], q)
                                      + ssm_smooth_eta_get_mean(eta[i], q));
      }
    }
    return draws;
}
