///
/// # Filtering Utility Functions
///
/// Functions used in filtering and log-likelihood calculations.
///
/// Requires: `ssm_utils.stan`
///

/**
Update the expected value of the predicted state, $a_{t + 1} = \E(\alpha_{t + 1} | \vec{y}_{1:t}$,

@param vector a An $m \times 1$ vector with the prected state, $a_t$.
@param vector c An $m \times 1$ vector with the system intercept, $c_t$
@param matrix T An $m \times m$ matrix with the transition matrix, $T_t$.
@param vector v A $p \times 1$ vector with the forecast error, $v_t$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $K_t$.
@return vector A $m \times 1$ vector with the predicted state at $t + 1$, $a_{t + 1}$.

The predicted state $a_{t + 1}$ is,
$$
a_{t + 1} = T_t a_t + K_t v_t + c_t .
$$

*/
vector ssm_filter_update_a(vector a, vector c, matrix T, vector v, matrix K) {
  vector[num_elements(a)] a_new;
  a_new = T * a + K * v + c;
  return a_new;
}

/**
Update the expected value of the predicted state, $P_{t + 1} = \Var(\alpha_{t + 1} | \vec{y}_{1:t}$,

@param matrix P An $m \times m$ vector with the variance of the prected state, $P_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $Z_t$.
@param matrix T An $m \times m$ matrix with the transition matrix, $T_t$.
@param matrix RQR A $m \times m$ matrix with the system covariance matrix, $R_t Q_t R_t'$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $K_t$.
@return matrix An $m \times 1$ vector with the predicted state at $t + 1$, $a_{t + 1}$.

The predicted state variance $P_{t + 1}$ is,
$$
P_{t + 1} = T_t P_t (T_t - K_t Z_t)' + R_t Q_t R_t' .
$$

*/
matrix ssm_filter_update_P(matrix P, matrix Z, matrix T,
                           matrix RQR, matrix K) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(T * P * (T - K * Z)' + RQR);
  return P_new;
}

/**
Update the forcast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}})$

@param matrix P An $m \times m$ vector with the variance of the prected state, $P_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $Z_t$.
@param matrix T An $m \times m$ matrix with the transition matrix, $T_t$.
@param matrix RQR An $m \times m$ matrix with the system covariance matrix, $R_t Q_t R_t'$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $K_t$.
@return vector An $m \times 1$ vector with the predicted state at $t + 1$, $a_{t + 1}$.

The forecast error $v_t$ is
$$
\vec{v}_t =\vec{y}_t - \mat{Z}_t \vec{a}_t - \vec{d}_t .
$$

*/
vector ssm_filter_update_v(vector y, vector a, vector d, matrix Z) {
  vector[num_elements(y)] v;
  v = y - Z * a - d;
  return v;
}

/**
Update the variance of the forcast error, $\mat{F}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))$

@param matrix P An $m \times m$ vector with the variance of the prected state, $P_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $Z_t$.
@param matrix H A $p \times p$ matrix with the observation covariance matrix, $H_t$.
@return matrix A $p \times p$ vector with $\mat{F}_t$.

The variance of the forecast error $\mat{F}_t$ is
$$
\mat{F}_t = \mat{Z}_t \mat{P}_t \mat{Z}_t + \mat{H}_t .
$$

*/
matrix ssm_filter_update_F(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] F;
  F = quad_form(P, Z') + H;
  return F;
}

/**
Update the precision of the forcast error, $\mat{F}^{-1}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))^{-1}$

@param matrix P An $m \times m$ vector with the variance of the prected state, $P_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $Z_t$.
@param matrix H A $p \times p$ matrix with the observation covariance matrix, $H_t$.
@return matrix A $p \times p$ vector with $\mat{F}^{-1}_t$.

This is the inverse of $\mat{F}_t$.

*/
matrix ssm_filter_update_Finv(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] Finv;
  Finv = inverse(ssm_filter_update_F(P, Z, H));
  return Finv;
}

/**
Update the Kalman gain, $\mat{K}_t$.

@param matrix P An $m \times m$ vector with the variance of the prected state, $P_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $Z_t$.
@param matrix T An $m \times m$ matrix with the transition matrix, $T_t$.
@param matrix Finv A $p \times p$ matrix
@return matrix An $m \times p$ matrix with the Kalman gain, $K_t$.

The Kalman gain is
$$
\mat{K}_t = \mat{T}_t \mat{P}_t \mat{Z}_t' F^{-1}_t .
$$

*/
matrix ssm_filter_update_K(matrix P, matrix Z, matrix T, matrix Finv) {
  matrix[cols(Z), rows(Z)] K;
  K = T * P * Z' * Finv;
  return K;
}

/**
Update $L_t$

@param matrix Z A $p \times m$ matrix with the design matrix, $Z_t$
@param matrix T An $m \times m$ matrix with the transition matrix, $T_t$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $K_t$.
@return matrix An $m \times m$ matrix, $L_t$.

$$
\mat{L}_t = \mat{T}_t - \mat{K}_t \mat{Z}_t .
$$

*/
matrix ssm_filter_update_L(matrix Z, matrix T, matrix K) {
  matrix[rows(T), cols(T)] L;
  L = T - K * Z;
  return L;
}

/**
Calculate the log-likelihood for a period

@param vector v A $p \times 1$ matrix with the forecast error, $v_t$.
@param matrix Finv A $p \times p$ matrix with variance of the forecast error, $F^{-1}_t$.
@return real An $m \times m$ matrix, $L_t$.

The log-likehood of a single observation in a state-space model is
$$
\ell_t = - \frac{1}{2} p \log(2 \pi) - \frac{1}{2} \left(\log|\mat{F}_t| + \vec{v}_t' \mat{F}^{-1}_t \vec{v}_t  \right)
$$
*/
real ssm_filter_update_ll(vector v, matrix Finv) {
  real ll;
  int p;
  p = num_elements(v);
  // det(A^{-1}) = 1 / det(A) -> log det(A^{-1}) = - log det(A)
  ll = (- 0.5 *
        (p * log(2 * pi())
         - log_determinant(Finv)
         + quad_form(Finv, v)
       ));
  return ll;
}
