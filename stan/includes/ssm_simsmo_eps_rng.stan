// ssm_simsmo_alpha_rng
vector[] ssm_simsmo_eps_rng(vector[] eps,
                      vector[] d, matrix[] Z, matrix[] H,
                      vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                      vector a1, matrix P1) {
    vector[dims(Z)[2]] draws[size(eps)];
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
        draws[i] = (ssm_sim_get_eps(sims[i], m, p, q)
                    - ssm_smooth_eps_get_mean(epshat_plus[i], p)
                    + ssm_smooth_eps_get_mean(eps[i], p));
      }
    }
    return draws;
}
