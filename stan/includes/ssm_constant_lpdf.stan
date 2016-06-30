///
/// # Constan Kalman Filter
///
/// Requires: `ssm_utils.stan`, `ssm_fitler_utils.stan`
///

/**
 Check if two matrices are approximately equal

 @param matrix A An $m \times n$ matrix.
 @param matrix B An $m \times n$ matrix.
 @param real The relative tolerance for convergence.

 The matrices $A$ and $B$ are considered approximately equal if
 $$
 \max(A - B) / \max(A) < \epsilon,
 $$
 where $\epsilon$ is the tolerance.

 */
int ssm_check_matrix_equal(matrix A, matrix B, real tol) {
  real eps;
  eps = max(to_vector(A - B)) / max(to_vector(A));
  if (eps < tol) {
    return 1;
  } else {
    return 0;
  }
}

/**
Log-likelihood of a Time-Invariant Linear Gaussian State Space Model

@param vector[] y Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
@param vector d Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
@param matrix Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param vector c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@param vector a1 Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
@param matrix P1 Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
@return real The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.

Unlike `ssm_filter`, this function requires the system matrices (`d`, `Z`, `H`, `c`, `T`, `R`, `Q`)
to all be time invariant (constant).
When the state space model is time-invariant, then the Kalman recursion for $\mat{P}_t$ converges.
This function takes advantage of this feature and stops updating $\mat{P}_t$ after it converges
to a steady state.

*/
real ssm_constant_lpdf(vector[] y,
                      vector d, matrix Z, matrix H,
                      vector c, matrix T, matrix R, matrix Q,
                      vector a1, matrix P1) {
  real ll;
  int n;
  int m;
  int p;

  n = size(y); // number of obs
  m = cols(Z);
  p = rows(Z);
  {
    vector[n] ll_obs;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    matrix[m, m] RQR;
    // indicator for if the filter has converged
    // This only works for time-invariant state space models
    int converged;
    matrix[m, m] P_old;
    real tol;
    converged = 0;
    tol = 1e-7;

    RQR = quad_form(Q, R);
    a = a1;
    P = P1;
    for (t in 1:n) {
      v = ssm_filter_update_v(y[t], a, d, Z);
      if (converged < 1) {
        Finv = ssm_filter_update_Finv(P, Z, H);
        K = ssm_filter_update_K(P, Z, T, Finv);
      }
      ll_obs[t] = ssm_filter_update_ll(v, Finv);
      // don't save a, P for last iteration
      if (t < n) {
        a = ssm_filter_update_a(a, c, T, v, K);
        // check for convergence
        // should only check for convergence if there are no missing values
        if (converged < 1) {
          P_old = P;
          P = ssm_filter_update_P(P, Z, T, RQR, K);
          converged = ssm_check_matrix_equal(P, P_old, tol);
        }
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}
