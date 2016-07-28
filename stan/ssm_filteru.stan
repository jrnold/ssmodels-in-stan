/**
`r render_section("Univariate Filtering")`

*/
/**
---
function: ssm_ufilter_idx
args:
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $6 \times 3$ integer array containing the indexes of the return values of the Kalman filter.
---

Indexes of the return values of the Kalman univariate filter function `ssm_ufilter`.

`ssm_ufilter_idx` returns a $6 \times 3$ integer array with the
(length, start index, stop index) of ($\ell_t$, $\vec{v}$, $\vec{F}^-1$, $\mat{K}$, $\vec{a}$, $\mat{P}$).

value            length                 start                               stop
---------------- ---------------------- ----------------------------------- ----------------------------------------------
$\ell_t$         $p$                    $1$                                 $p$
$\vec{v}$        $p$                    $1 + p$                             $2 p$
$\mat{F}^{-1}$   $p$                    $1 + 2 p$                           $3 p$
$\mat{K}$        $mp$                   $1 + 3 p$                           $3 p + mp$
$\vec{a}_t$      $m$                    $1 + 3 p + mp$                      $3 p + mp + m$
$\mat{P}^t$      $m (m + 1) / 2$        $1 + 3 p + mp + m$                  $3 p + mp + m (m + 1) / 2$


*/

int[,] ssm_ufilter_idx(int m, int p) {
  int sz[6, 3];
  // loglike
  sz[1, 1] = p;
  // v
  sz[2, 1] = p;
  // Finv
  sz[3, 1] = p;
  // K
  sz[4, 1] = m * p;
  // a
  sz[5, 1] = m;
  // P
  sz[6, 1] = symmat_size(m);
  // Fill in start and stop points
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:6) {
    sz[i, 2] = sz[i - 1, 3] + 1;
    sz[i, 3] = sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}

/**
---
function: ssm_ufilter_size
args:
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: The number of elements in the vector.
---

Number of elements in vector containing the results from `ssm_ufilter`.


*/

int ssm_ufilter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  sz = idx[6, 3];
  return sz;
}

/**
---
function: ssm_ufilter_get_loglik
args:
- name: x
  description: vector with results from `ssm_ufilter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: Vector of length $p$ with the log-likelihood, $\vec{\ell}_t$, of each observation.
---

Get the log-likehood from the results of `ssm_ufilter`.


*/

vector ssm_ufilter_get_loglik(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  y = segment(x, idx[1, 2], idx[1, 1]);
  return y;
}

/**
---
function: ssm_ufilter_get_v
args:
- name: x
  description: vector with results from `ssm_ufilter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $p \times 1$ vector with the forecast error, $\vec{v}_t$.
---

Get the forecast error from the results of `ssm_ufilter`.


*/

vector ssm_ufilter_get_v(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  y = segment(x, idx[2, 2], idx[2, 1]);
  return y;
}

/**
---
function: ssm_ufilter_get_Finv
args:
- name: x
  description: vector with results from `ssm_ufilter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A length $p$ vector with the forecast precisions, $\vec{f}^{-1}_t$.
---

Get the forecast precision from the results of `ssm_ufilter`.


*/

vector ssm_ufilter_get_Finv(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  y = segment(x, idx[3, 2], idx[3, 1]);
  return y;
}

/**
---
function: ssm_ufilter_get_K
args:
- name: x
  description: vector with results from `ssm_ufilter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
---

Get the Kalman gain from the results of `ssm_ufilter`.


*/

matrix ssm_ufilter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[4, 2], idx[4, 1]), m, p);
  return y;
}

/**
---
function: ssm_ufilter_get_a
args:
- name: A
  description: vector with results from `ssm_ufilter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: An $m \times 1$ vector with the expected value of the predicted state, $\E(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \vec{a}_t$.
---

Get the expected value of the predicted state from the results of `ssm_ufilter`.


*/

vector ssm_ufilter_get_a(vector x, int m, int p) {
  vector[m] y;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  y = segment(x, idx[5, 2], idx[5, 1]);
  return y;
}

/**
---
function: ssm_ufilter_get_P
args:
- name: A
  description: vector with results from `ssm_ufilter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: An $m \times m$ matrix with the variance of the predicted state, $\Var(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \mat{P}_t$.
---

Get the variance of the predicted state from the results of `ssm_ufilter`.


*/

matrix ssm_ufilter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  int idx[6, 3];
  idx = ssm_ufilter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[6, 2], idx[6, 1]), m);
  return y;
}

/**
---
function: ssm_ufilter
args:
- name: y
  description: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Diagonal of the observation covariance matrix, $\mat{H}_t = \diag(\vec{h_t})$. An array of $p \times 1$ vectors.
- name: c
  description: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- name: T
  description: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- name: R
  description: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- name: Q
  description: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- name: a1
  description: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- name: P1
  description: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
returns: Array of size $n$ of $(3 p + mp + m + m (m + 1) / 2) \times 1$ vectors in the format described in `ssm_ufilter_idx`.
---


Kalman filter


For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

`ssm_ufilter` runs a forward filter on the state space model and calculates,

- log-likelihood for each observation, $\ell_t$.
- Forecast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y}_{1:(t -1)})$.
- Forecast precision, $\mat{F}^{-1}_t$.
- Kalman gain, $\mat{K}_t$.
- Predicted states, $\vec{a}_t = \E(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.
- Variance of the predicted states, $\mat{P}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.

The results of Kalman filter for a given are returned as a $1 + p + p (p + 1) / 2 + m p + m (m + 1) / 2$ vector for each time period, where
$$
(\ell_t, \vec{v}_t', \VEC(\mat{F}^{-1}_t)', \VEC(\mat{K}_t)', \vec{a}_t', \VEC(\mat{P}_t)' )'.
$$

*/

vector[] ssm_ufilter(vector[] y,
                    vector[] d, matrix[] Z, vector[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {

  // returned data
  vector[ssm_ufilter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
  int q;
  int n;
  int p;
  int m;

  // sizes
  n = size(y); // number of obs
  m = dims(Z)[3]; // number of states
  p = dims(Z)[2]; // obs size
  q = dims(Q)[2]; // number of state disturbances

  //print("Sizes: n = ", m, ", p = ", n, ", m = ", m, ", q = ", q);
  {
    // system matrices for current iteration
    vector[p] d_t;
    real d_ti;
    matrix[p, m] Z_t;
    row_vector[m] Z_ti;
    vector[p] H_t;
    real h_ti;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    // result matricees for each iteration
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    vector[p] Finv;
    vector[m] K_i;
    matrix[m, p] K;
    vector[p] ll;
    int idx[6, 3];

    idx = ssm_ufilter_idx(m, p);

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form_sym(Q_t, R_t ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      # Save predictes for a_{t,1} and P_{t,1}
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);

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
        if (size(R) > 1 || size(Q) > 1) {
          RQR = quad_form_sym(Q_t, R_t ');
        }
      }
      // updating
      for (i in 1:p) {
        Z_ti = row(Z_t, i);
        v[i] = ssm_update_v_u(y[t, i], a, d_t[i], Z_ti);
        Finv[i] = ssm_update_Finv_u(P, Z_ti, H_t[i]);
        K_i = ssm_update_K_u(P, Z_ti, Finv);
        K[:, i] = K_i;
        ll[i] = ssm_update_loglik_u(v[i], Finv[i]);
        a = ssm_update_a_u1(a, v[i], K_i);
        P = ssm_update_P_u1(P, Finv[i], K_i);
      }
      // saving
      res[t, idx[1, 2]:idx[1, d]]] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = Finv;
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a = ssm_update_a_u2(a, c_t, T_t);
        P = ssm_update_P_u2(P, T_t, RQR);
      }
    }
  }
  return res;
}


vector[] ssm_ufilter_miss(vector[] y,
                    vector[] d, matrix[] Z, vector[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1, int[,] miss) {

  // returned data
  vector[ssm_ufilter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
  int q;
  int n;
  int p;
  int m;

  // sizes
  n = size(y); // number of obs
  m = dims(Z)[3]; // number of states
  p = dims(Z)[2]; // obs size
  q = dims(Q)[2]; // number of state disturbances

  //print("Sizes: n = ", m, ", p = ", n, ", m = ", m, ", q = ", q);
  {
    // system matrices for current iteration
    vector[p] d_t;
    real d_ti;
    matrix[p, m] Z_t;
    row_vector[m] Z_ti;
    vector[p] H_t;
    real h_ti;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    // result matricees for each iteration
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    vector[p] Finv;
    vector[m] K_i;
    matrix[m, p] K;
    vector[p] ll;
    int idx[6, 3];

    idx = ssm_ufilter_idx(m, p);

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form_sym(Q_t, R_t ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      # Save predictes for a_{t,1} and P_{t,1}
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);

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
        if (size(R) > 1 || size(Q) > 1) {
          RQR = quad_form_sym(Q_t, R_t ');
        }
      }
      // updating
      for (i in 1:p) {
        if (miss[t, i]) {
          v[i] = 0.;
          Finv[i] = 0.;
          K[:, i] = rep_vector(0., m);
          ll[i] = 0.;
        } else {
          Z_ti = row(Z_t, i);
          v[i] = ssm_update_v_u(y[t, i], a, d_t[i], Z_ti);
          Finv[i] = ssm_update_Finv_u(P, Z_ti, H_t[i]);
          K_i = ssm_update_K_u(P, Z_ti, Finv);
          K[:, i] = K_i;
          ll[i] = ssm_update_loglik_u(v[i], Finv[i]);
          a = ssm_update_a_u1(a, v[i], K_i);
          P = ssm_update_P_u1(P, Finv[i], K_i);
        }
      }
      // saving
      res[t, idx[1, 2]:idx[1, d]]] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = Finv;
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a = ssm_update_a_u2(a, c_t, T_t);
        P = ssm_update_P_u2(P, T_t, RQR);
      }
    }
  }
  return res;
}

/**
---
function: ssm_ufilter_states
args:
- name: m
  description: Number of states
returns: The size of the vector
---

Length of the vectors returned by `ssm_ufilter_states`


*/

int ssm_ufilter_states_size(int m) {
  return ssm_filter_states_size(m);
}

/**
---
function: ssm_ufilter_states_get_a
args:
- name: x
  description: A vector returned by `ssm_ufilter_states`
- name: m
  description: Number of states
returns: An $m \times 1$ vector with the filtered expected value of the state, $\vec{a}_{t|t} = \E(\vec{\alpha}_t | \vec{y}_{1:t})$.
---

Extract $a_{t|t}$ from the results of `ssm_ufilter_states`


*/

vector ssm_ufilter_states_get_a(vector x, int m) {
  return ssm_filter_states_get_a(x, m);
}

/**
---
function: ssm_ufilter_states_get_P
args:
- name: x
  description: A vector returned by `ssm_ufilter_states`
- name: m
  description: Number of states
returns: An $m \times m$ matrix with the filtered variance of the state, $\mat{P}_{t|t} = \Var(\vec{\alpha}_t | \vec{y}_{1:t})$.
---

Extract $P_{t|t}$ from the results of `ssm_ufilter_states`


*/

matrix ssm_ufilter_states_get_P(vector x, int m) {
  return ssm_filter_states_get_P(x, m);
}

/**
---
function: ssm_ufilter_states_update_a
args:
- name: a
  description: An $m \times 1$ vector with the expected value of the predicted state, $\vec{a}_t$.
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: v
  description: A $p \times 1$ vector with the forecast errors, $\vec{v}_t$
- name: Finv
  description: A $p \times 1$ vector with the forecast prediction, $\vec{f}_{t}^{-1}$.
returns: An $m \times 1$ matrix the expected value of the fitered states, $\E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t}$.
---


Calculate filtered state values [@DurbinKoopman2012, Sec 4.3.2], from the results of
a univariate filter (`ssm_`ufilter` or `ssm_ufilter_miss`):
$$
\E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t} = ....
$$


*/

vector ssm_ufilter_states_update_a(vector a, matrix P, matrix Z,
                                  vector v, vector Finv) {
  vector[num_elements(a)] aa;
  aa = a + P * Z ' * Finv * v;
  return aa;
}

/**
---
function: ssm_ufilter_states_update_P
args:
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: Finv
  description: A $p \times p$ matrix with the forecast prediction, $\mat{F}_{t}^{-1}$.
returns: An $m \times m$ matrix with variance of the filtered states, $\Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t}$.
---


Calculate filtered state variance values [@DurbinKoopman2012, Sec 4.3.2], from the reults of
a univariate filter (`ssm_ufilter` or `ssm_ufilter_miss`):
$$
\Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t} = ...
$$



*/

matrix ssm_ufilter_states_update_P(matrix P, matrix Z, vector Finv) {
  matrix[rows(P), cols(P)] PP;
  PP = to_symmetric_matrix(P - P * quad_form_sym(Finv, Z) * P);
  return PP;
}


/**
---
function: ssm_ufilter_states
args:
- name: filter
  description: Results from `ssm_ufilter`
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
returns: of size $n$ of vectors.
---

Calculate filtered expected values and variances of the states

The filtering function `ssm_ufilter` returns the mean and variance of the predicted states,
$\vec{a}_t = \E(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$ and $\mat{P}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.


The vectors returned by `ssm_ufilter_states` are of length $m + m ^ 2$, with
$$
\vec{v}_t = (\vec{a}_{t|t}', \VEC(\vec{P}_{t|t})' )'
$$
Use the functions `ssm_ufilter_states_get_a` and `ssm_ufilter_states_get_P` to extract
elements from the results.

For `Z` the array can have a size of 1, if it is not time-varying, or a size of $n - 1$ if it is time varying.

*/

vector[] ssm_ufilter_states(vector[] filter, matrix[] Z) {
  vector[ssm_ufilter_states_size(dims(Z)[3])] res[size(filter)];
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
    vector[p] v;
    vector[p] Finv;
    matrix[m, p] K;
    vector[m] K_i;
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
      v = ssm_ufilter_get_v(filter[t], m, p);
      Finv = ssm_ufilter_get_Finv(filter[t], m, p);
      K = ssm_ufilter_get_K(filter[t], m, p);
      a = ssm_ufilter_get_a(filter[t], m, p);
      P = ssm_ufilter_get_P(filter[t], m, p);
      for (i in 1:p) {
        K_i = col(K, i);
        a = ssm_update_a_u1(a, v[i], K_i);
        P = ssm_update_P_u1(P, Finv[i], K_i);
      }
      res[t, :m] = a;
      res[t, (m + 1): ] = symmat_to_vector(P);
    }
  }
  return res;
}
