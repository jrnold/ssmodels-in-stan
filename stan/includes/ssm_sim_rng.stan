///
/// # Simulators and Smoothing Simulators
///

/**
Indexes of each component of `ssm_sim_rng` results.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return A 4 x 3 array of integers with the (length, start location, and end location)
  of $y_t$, $\alpha_t$, $\varepsilon_t$, and $\eta_t$ in the results of `ssm_sim_rng`.

element         length         start         end
--------------- -------------- ------------- -----------
$y_t$           $p$            $1$           $p$
$\alpha$_t      $m$            $p + 1$       $p + m$
$\varepsilon_t$ $p$            $p + m + 1$   $2 p + m$
$\eta_t$        $q$            $2 p + m + 1$ $2 p + m + q$

It is preferrable to use `ssm_sim_get_y`, `ssm_sim_get_a`, `ssm_sim_get_eps`,
and `ssm_sim_get_eta` to extract values from these vectors.

*/
int[,] ssm_sim_idx(int m, int p, int q) {
  int sz[4, 3];
  // y
  sz[1, 1] = p;
  // a
  sz[2, 1] = m;
  // eps
  sz[3, 1] = p;
  // eta
  sz[4, 1] = q;
  // Fill in start and stop points
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:4) {
    sz[i, 2] = sz[i - 1, 3] + 1;
    sz[i, 3] = sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}

/**
The number of elements in vectors returned by `ssm_sim_rng` results.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return int The number of elements

*/
int ssm_sim_size(int m, int p, int q) {
  int sz;
  sz = ssm_sim_idx(m, p, q)[4, 3];
  return sz;
}

/**
Extract $\vec{y}_t$ from vectors returne by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector vector A $p \times 1$ vector with $\vec{y}_t$.

*/
vector ssm_sim_get_y(vector x, int m, int p, int q) {
  vector[m] y;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  y = x[idx[1, 2]:idx[1, 3]];
  return y;
}

/**
Extract $\vec{\alpha}_t$ from vectors returne by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector A $m \times 1$ vector with $\vec{\alpha}_t$.

*/
vector ssm_sim_get_a(vector x, int m, int p, int q) {
  vector[m] a;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  a = x[idx[2, 2]:idx[2, 3]];
  return a;
}

/**
Extract $\vec{\varepsilon}_t$ from vectors returne by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector vector A $p \times 1$ vector with $\vec{\varepsilon}_t$.


*/
vector ssm_sim_get_eps(vector x, int m, int p, int q) {
  vector[m] eps;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  eps = x[idx[3, 2]:idx[3, 3]];
  return eps;
}

/**
Extract $\vec{\eta}_t$ from vectors returne by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector vector A $q \times 1$ vector with $\vec{\eta}_t$.

*/
vector ssm_sim_get_eta(vector x, int m, int p, int q) {
  vector[m] eta;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  eta = x[idx[4, 2]:idx[4, 3]];
  return eta;
}

/**
Simulate from a Linear Gaussian State Space model.

@param vector[] y Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
@param vector[] d Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param vector[] c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@param vector a1 Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
@param matrix P1 Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
@return Array of size $n$ of vectors with Draw $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\eta}_t$ and $\vec{\varepsilon}}_t$. See the description.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

Draw $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\eta}_t$ and $\vec{\varepsilon}}_t$ from
the state space model,
$$
\begin{aligned}[t]
\vec{y}_t &= \vec{d}_t + \mat{Z}_t \vec{\alpha}_t + \vec{\varepsilon}_t,  &
\vec{\varepsilon}_t & \sim N(0, \mat{H}_t), \\
\vec{\alpha}_{t + 1} &= \vec{c}_t + \mat{T}_t \vec{\alpha}_t + \mat{R}_t \vec{\eta}_t,  &
\vec{\eta}_t & \sim N(0, \mat{Q}_t), \\
&& \vec{\alpha}_1 &\sim N(\vec{a}_1, \mat{P}_1) .
\end{aligned}
$$

The returned vectors are of length $2 p + m + q$, in the format,
$$
(\vec{y}_t', \vec{\alpha}_t', \vec{\varepsilon}_t', \vec{\eta}_t') .
$$
Note that $\eta_n = \vec{0}_q$.
Use the functions `ssm_sim_get_y`, `ssm_sim_get_a`, `ssm_sim_get_eps`, and
`ssm_sim_get_eta` to extract values from the vector.


*/
vector[] ssm_sim_rng(int n,
                    vector[] d, matrix[] Z, matrix[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {
  vector[ssm_sim_size(dims(Z)[3], dims(Z)[2], dims(Q)[2])] ret[n];
  int p;
  int m;
  int q;
  p = dims(Z)[2];
  m = dims(Z)[3];
  q = dims(Q)[2];
  {
    // system matrices for current iteration
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    // outputs
    vector[p] y;
    vector[p] eps;
    vector[m] a;
    vector[q] eta;
    // constants
    vector[p] zero_p;
    vector[q] zero_q;
    vector[m] zero_m;
    int idx[4, 3];

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];

    idx = ssm_sim_idx(m, p, q);
    zero_p = rep_vector(0.0, p);
    zero_q = rep_vector(0.0, q);
    zero_m = rep_vector(0.0, m);
    a = multi_normal_rng(a1, P1);
    for (t in 1:n) {
      // set system matrices
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
        }
        // system matrices are n - 1 length
        if (t < n) {
          if (size(c) > 1) {
            c_t = c[t];
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
        }
      }
      // draw forecast error
      eps = multi_normal_rng(zero_p, H_t);
      // draw observed value
      y = d_t + Z_t * a + eps;
      // since eta_t is for alpha_{t + 1}, we don't
      // draw it for t == n
      if (t == n) {
        eta = zero_q;
      } else {
        eta = multi_normal_rng(zero_q, Q_t);
      }
      // save
      ret[t, idx[1, 2]:idx[1, 3]] = y;
      ret[t, idx[2, 2]:idx[2, 3]] = a;
      ret[t, idx[3, 2]:idx[3, 3]] = eps;
      ret[t, idx[4, 2]:idx[4, 3]] = eta;
      // a_{t + 1}
      if (t < n) {
        a = c_t + T_t * a + R_t * eta;
      }
    }
  }
  return ret;
}
