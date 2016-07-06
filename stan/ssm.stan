/**
@chapter Stan Functions

State space functionality for Stan is provided as a set of user-defined functions.

Add the following line to the Stan model file in which depends on these functions.
```stan
functions {
  #include ssm.stan
  // other functions ...
}
```

To actually include the functions in the model, you need to use the function `stanc_builder`,
instead of `stan` or `stanc`:
```r
model <- stanc_builder("yourmodel.stan", isystem = "path/to/ssm/")
stan(model_code = model$model_code)
```


@section Utility Functions

*/
/** to_symmetric_matrix
Ensure a matrix is symmetrix.

@param x An $n \times n$ matrix.
@return An $n \times n$ symmetric matrix: $0.5 (x + x')$.

*/
matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}

/** to_matrix_colwise

Convert vector to a matrix (column-major).

@param vector v An $n \times m$ vector.
@param int m Number of rows in the vector
@param int n Number of columns in the vector
@return matrix A $m \times n$ matrix containting the elements from `v`

*/
matrix to_matrix_colwise(vector v, int m, int n) {
  matrix[m, n] res;
  int k;
  k = 1;
  // by col
  for (j in 1:n) {
    // by row
    for (i in 1:m) {
      res[i, j] = v[k];
      k = k + 1;
    }
  }
  return res;
}

/** symmat_size
Calculate the number of unique elements in a symmetric matrix

The number of unique elements in an $m \times m$ matrix is
$(m \times (m + 1)) / 2$.

@param matrix x An $m \times m$ matrix.
@return int The number of unique elements

*/
int symmat_size(int n) {
  int sz;
  // This calculates it iteratively because Stan gives a warning
  // with integer division.
  sz = 0;
  for (i in 1:n) {
    sz = sz + i;
  }
  return sz;
}

/** find_symmat_dim

Given vector with $n$ elements containing the $m (m + 1) / 2$ elements of a symmetric matrix,
return $m$.

@param int n The number of unique elements in a symmetric matrix.
@return int The dimension of the associated symmetric matrix.

*/
int find_symmat_dim(int n) {
  // This could be solved by finding the positive root of $m = m (m + 1)/2 but
  // Stan doesn't support all the functions necessary to do this.
  int i;
  int remainder;
  remainder = n;
  i = 0;
  while (remainder > 0) {
    i = i + 1;
    remainder = remainder - i;
  }
  return i;
}

/** vector_to_symmat
Convert a vector to a symmetric matrix

@param vector x The vector with the unique elements
@param int n The dimensions of the returned matrix: $n \times n$.
@return matrix An $n \times n$ symmetric matrix.

*/
matrix vector_to_symmat(vector x, int n) {
  matrix[n, n] m;
  int k;
  k = 1;
  // for column
  for (j in 1:n) {
    // for row
    for (i in j:n) {
      m[i, j] = x[k];
      if (i != j) {
        m[j, i] = m[i, j];
      }
      k = k + 1;
    }
  }
  return m;
}

/** symmat_to_vector
Convert an $n \times n$ symmetric matrix to a length $n (n + 1) / 2$ vector
containing its unique elements.

@param vector x An $n \times n$ matrix.
@return vector A $n (n + 1) / 2$ vector with the unique elements in $x$.

*/
vector symmat_to_vector(matrix x) {
  vector[symmat_size(min(rows(x), cols(x)))] v;
  int m;
  int k;
  k = 1;
  m = min(rows(x), cols(x));
  // if x is m x n symmetric, then this will return
  // only parts of an m x m matrix.
  for (j in 1:m) {
    for (i in j:m) {
      v[k] = x[i, j];
      k = k + 1;
    }
  }
  return v;
}
/**

@section Filtering

Functions used in filtering and log-likelihood calculations.
*/
/** ssm_filter_update_a
Update the expected value of the predicted state, $\vec{a}_{t + 1} = \E(\vec{\alpha}_{t + 1} | \vec{y}_{1:t})$,

@param vector a An $m \times 1$ vector with the predicted state, $\vec{a}_t$.
@param vector c An $m \times 1$ vector with the system intercept, $\vec{c}_t$
@param matrix T An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
@param vector v A $p \times 1$ vector with the forecast error, $\vec{v}_t$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
@return vector A $m \times 1$ vector with the predicted state at $t + 1$, $\vec{a}_{t + 1}$.

The predicted state $\vec{a}_{t + 1}$ is,
$$
\vec{a}_{t + 1} = \mat{T}_t \vec{a}_t + \mat{K}_t \vec{v}_t + \vec{c}_t .
$$

*/
vector ssm_filter_update_a(vector a, vector c, matrix T, vector v, matrix K) {
  vector[num_elements(a)] a_new;
  a_new = T * a + K * v + c;
  return a_new;
}

/** ssm_filter_update_P
Update the expected value of the predicted state, $\mat{P}_{t + 1} = \Var(\alpha_{t + 1} | \vec{y}_{1:t})$,

@param matrix P An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param matrix T An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
@param matrix RQR A $m \times m$ matrix with the system covariance matrix, $\mat{R}_t \mat{Q}_t \mat{R}_t'$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
@return matrix An $m \times 1$ vector with the predicted state at $t + 1$, $\vec{a}_{t + 1}$.

The predicted state variance $\mat{P}_{t + 1}$ is,
$$
\mat{P}_{t + 1} = \mat{T}_t \mat{P}_t (\mat{T}_t - \mat{K}_t \mat{Z}_t)' + \mat{R}_t \mat{Q}_t \mat{R}_t' .
$$

*/
matrix ssm_filter_update_P(matrix P, matrix Z, matrix T,
                           matrix RQR, matrix K) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(T * P * (T - K * Z)' + RQR);
  return P_new;
}

/** ssm_filter_update_v
Update the forcast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}})$

@param matrix P An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param matrix T An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
@param matrix RQR An $m \times m$ matrix with the system covariance matrix, $\mat{R}_t \mat{Q}_t \mat{R}_t'$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
@return vector An $m \times 1$ vector with the predicted state at $t + 1$, $\vec{a}_{t + 1}$.

The forecast error $\vec{v}_t$ is
$$
\vec{v}_t =\vec{y}_t - \mat{Z}_t \vec{a}_t - \vec{d}_t .
$$

*/
vector ssm_filter_update_v(vector y, vector a, vector d, matrix Z) {
  vector[num_elements(y)] v;
  v = y - Z * a - d;
  return v;
}

/** ssm_filter_update_F
Update the variance of the forcast error, $\mat{F}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))$

@param matrix P An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param matrix H A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.
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

/** ssm_filter_update_Finv
Update the precision of the forcast error, $\mat{F}^{-1}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))^{-1}$

@param matrix P An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param matrix H A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.
@return matrix A $p \times p$ vector with $\mat{F}^{-1}_t$.

This is the inverse of $\mat{F}_t$.

*/
matrix ssm_filter_update_Finv(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] Finv;
  Finv = inverse(ssm_filter_update_F(P, Z, H));
  return Finv;
}

/** ssm_filter_update_K
Update the Kalman gain, $\mat{K}_t$.

@param matrix P An $m \times m$ vector with the variance of the predicted state, $P_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param matrix T An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
@param matrix Finv A $p \times p$ matrix
@return matrix An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.

The Kalman gain is
$$
\mat{K}_t = \mat{T}_t \mat{P}_t \mat{Z}_t' \mat{F}^{-1}_t .
$$

*/
matrix ssm_filter_update_K(matrix P, matrix Z, matrix T, matrix Finv) {
  matrix[cols(Z), rows(Z)] K;
  K = T * P * Z' * Finv;
  return K;
}

/** ssm_filter_update_L
Update $L_t$

@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$
@param matrix T An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
@param matrix K An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
@return matrix An $m \times m$ matrix, $\mat{L}_t$.

$$
\mat{L}_t = \mat{T}_t - \mat{K}_t \mat{Z}_t .
$$

*/
matrix ssm_filter_update_L(matrix Z, matrix T, matrix K) {
  matrix[rows(T), cols(T)] L;
  L = T - K * Z;
  return L;
}

/** ssm_filter_update_ll
Calculate the log-likelihood for a period

@param vector v A $p \times 1$ matrix with the forecast error, $\vec{v}_t$.
@param matrix Finv A $p \times p$ matrix with variance of the forecast error, $\mat{F}^{-1}_t$.
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

/**
@section Filtering

*/
/** ssm_filter_idx
Indexes of the return values of the Kalman filter functions:
`ssm_filter`.

@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return int[,] A $6 \times 3$ integer array containing the indexes of the return values of the Kalman filter.

`ssm_filter_idx` returns a $6 \times 3$ integer array with the
(length, start index, stop index) of ($\ell_t$, $\vec{v}$, $\vec{F}^-1$, $\mat{K}$, $\vec{a}$, $\mat{P}$).

value            length                 start                               stop
---------------- ---------------------- ----------------------------------- ----------------------------------------------
$\ell_t$         $1$                    $1$                                 $1$
$\vec{v}$        $p$                    $2$                                 $1 + p$
$\mat{F}^{-1}$   $p (p + 1) / 2$        $2 + p$                             $1 + p + p * (p + 1) / 2$
$\mat{K}$        $mp$                   $2 + p + p (p + 1) / 2$             $1 + p + p (p + 1) / 2 + mp$
$\vec{a}_t$      $m$                    $2 + p + p (p + 1) / 2 + mp$        $1 + p + p (p + 1) / 2 + mp + m$
$\mat{P}^t$      $m (m + 1) / 2$        $2 + p + p (p + 1) / 2 + mp + m$    $1 + p + p (p + 1) / 2 + mp + m * (m + 1) / 2$


*/
int[,] ssm_filter_idx(int m, int p) {
  int sz[6, 3];
  // loglike
  sz[1, 1] = 1;
  // v
  sz[2, 1] = p;
  // Finv
  sz[3, 1] = symmat_size(p);
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

/** ssm_filter_size
Number of elements in vector containing filter results

@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return int The number of elements in the vector.

*/
int ssm_filter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  sz = idx[6, 3];
  return sz;
}

/** ssm_filter_get_loglik
Get the log-likehood from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return real The log-likelihood $\ell_t$

*/
real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y = x[1];
  return y;
}

/** ssm_filter_get_v
Get the forecast error from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return vector A $p \times 1$ vector with the forecast error, $\vec{v}_t$.

*/
vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[2, 2], idx[2, 1]);
  return y;
}

/** ssm_filter_get_Finv
Get the forecast precision from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return matrix A $p \times p$ matrix with the forecast precision, $\mat{F}^{-1}_t$.

*/
matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[3, 2], idx[3, 1]), p);
  return y;
}

/** ssm_filter_get_K
Get the Kalman gain from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return matrix A $m \times p$ matrix with the Kalman gain, $\mat{F}^{-1}_t$.

*/
matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[4, 2], idx[4, 1]), m, p);
  return y;
}

/** ssm_filter_get_a
Get the expected value of the predicted state from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return vector An $m \times 1$ vector with the expected value of the predicted state, $\E(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \vec{a}_t$.

*/
vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[5, 2], idx[5, 1]);
  return y;
}

/** ssm_filter_get_P
Get the variance of the predicted state from the results of `ssm_filter`.

@param vector A vector with results from `ssm_filter`.
@param int m The number of states
@param int p The size of the observation vector $\vec{y}_t$.
@return matrix An $m \times m$ matrix with the variance of the predicted state, $\Var(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \mat{P}_t$.

*/
matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[6, 2], idx[6, 1]), m);
  return y;
}

/** ssm_filter

Kalman filter

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
@return vector[] Array of size $n$ of $(1 + p + p (p + 1) / 2 + mp + m + m (m + 1) / 2) \times 1$ vectors in the format described in `ssm_filter_idx`.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

`ssm_filter` runs a forward filter on the state space model and calculates,

- log-likelihood for each observation, $\ell_t$.
- Forecast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y}_{1:(t -1)})$.
- Forecast precision, $\mat{F}^{-1}_t$.
- Kalman gain, $\mat{K}_t$.
- Predicted states, $\vec{a}_t = \E(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.
- Variance of the predicted states, $\mat{P}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.

The results of Kalman filter for a given are returned as a $1 + p + p * (p + 1) / 2 + m p + m * (m + 1) / 2$ vector for each time period, where
$$
(\ell_t, \vec{v}_t', \VEC(\mat{F}^{-1}_t)', \VEC(\mat{K}_t)', \vec{a}_t', \VEC(\mat{P}_t)' )'.
$$

*/
vector[] ssm_filter(vector[] y,
                    vector[] d, matrix[] Z, matrix[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {

  // returned data
  vector[ssm_filter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
  int q;
  int n;
  int p;
  int m;

  // sizes
  n = size(y); // number of obs
  p = dims(Z)[2]; // obs size
  m = dims(Z)[3]; // number of states
  q = dims(Q)[2]; // number of state disturbances

  //print("Sizes: n = ", m, ", p = ", n, ", m = ", m, ", q = ", q);
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
    // result matricees for each iteration
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    real ll;
    int idx[6, 3];

    idx = ssm_filter_idx(m, p);

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form(Q_t, R_t);
    a = a1;
    P = P1;
    for (t in 1:n) {
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
        if (size(R) > 1 && size(Q) > 1) {
          RQR = quad_form(Q_t, R_t);
        }
      }
      // updating
      v = ssm_filter_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_filter_update_Finv(P, Z_t, H_t);
      K = ssm_filter_update_K(P, T_t, Z_t, Finv);
      ll = ssm_filter_update_ll(v, Finv);
      // saving
      res[t, 1] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = symmat_to_vector(Finv);
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a = ssm_filter_update_a(a, c_t, T_t, v, K);
        P = ssm_filter_update_P(P, Z_t, T_t, RQR, K);
      }
    }
  }
  return res;
}

/** ssm_filter_states
Length of the vectors returned by `ssm_filter_states`

@param int m Number of states
@return int The size of the vector

*/
int ssm_filter_states_size(int m) {
  int sz;
  sz = m + symmat_size(m);
  return sz;
}

/** ssm_filter_states_get_a
Extract $a_{t|t}$ from the results of `ssm_filter_states`

@param vector x A vector returned by `ssm_filter_states`
@param int m Number of states
@return matrix An $m \times 1$ vector with the filtered expected value of the state, $\vec{a}_{t|t} = \E(\vec{\alpha}_t | \vec{y}_{1:t})$.

*/
vector ssm_filter_states_get_a(vector x, int m) {
  vector[m] a;
  a = x[ :m];
  return a;
}

/** ssm_filter_states_get_P
Extract $P_{t|t}$ from the results of `ssm_filter_states`

@param vector x A vector returned by `ssm_filter_states`
@param int m Number of states
@return matrix An $m \times m$ matrix with the filtered variance of the state, $\mat{P}_{t|t} = \Var(\vec{\alpha}_t | \vec{y}_{1:t})$.

*/
matrix ssm_filter_states_get_P(vector x, int m) {
  matrix[m, m] P;
  P = vector_to_symmat(x[(m + 1): ], m);
  return P;
}

/** ssm_filter_states_update_a

Calculate filtered state values [@DurbinKoopman2012, Sec 4.3.2],
$$
\E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t} = \mat{T}_t * \vec{a}_t + \mat{K}_t \vec{v}_t .
$$

@param vector a An $m \times 1$ vector with the expected value of the predicted state, $\vec{a}_t$.
@param matrix P An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param vector v A $p \times 1$ vector with the forecast errors, $\vec{v}_t$
@param matrix Finv A $p \times p$ matrix with the forecast prediction, $\mat{F}_{t}^{-1}$.
@return matrix An $m \times 1$ matrix the expected value of the fitered states, $\E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t}$.

*/
vector ssm_filter_states_update_a(vector a, matrix P, matrix Z,
                                  vector v, matrix Finv) {
  vector[num_elements(a)] aa;
  aa = a + P * Z ' * Finv * v;
  return aa;
}

/** ssm_filter_states_update_P

Calculate filtered state variance values [@DurbinKoopman2012, Sec 4.3.2],
$$
\Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t} = \mat{P}_t - \mat{P}_t \mat{Z}_t' \mat{F}_t^{-1} \mat{Z}_t \mat{P}_t .
$$

@param matrix P An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
@param matrix Z A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
@param matrix Finv A $p \times p$ matrix with the forecast prediction, $\mat{F}_{t}^{-1}$.
@return matrix An $m \times m$ matrix with variance of the filtered states, $\Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t}$.


*/
matrix ssm_filter_states_update_P(matrix P, matrix Z, matrix Finv) {
  matrix[rows(P), cols(P)] PP;
  PP = to_symmetric_matrix(P - P * quad_form(Finv, Z) * P);
  return PP;
}


/** ssm_filter_states
Calculate filtered expected values and variances of the states

The filtering function `ssm_filter` returns the mean and variance of the predicted states,
$\vec{a}_t = \E(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$ and $\mat{P}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.

@param vector[] filter Results from `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@return Array of size $n$ of vectors.

The vectors returned by `ssm_filter_states` are of length $m + m ^ 2$, with
$$
\vec{v}_t = (\vec{a}_{t|t}', \VEC(\vec{P}_{t|t})' )'
$$
Use the functions `ssm_filter_states_get_a` and `ssm_filter_states_get_P` to extract
elements from the results.

For `Z` the array can have a size of 1, if it is not time-varying, or a size of $n - 1$ if it is time varying.

*/
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
      // extract values from the filter
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      a = ssm_filter_get_a(filter[t], m, p);
      P = ssm_filter_get_P(filter[t], m, p);
      // calcualte filtered values
      aa = ssm_filter_states_update_a(a, P, Z_t, v, Finv);
      PP = ssm_filter_states_update_P(P, Z_t, Finv);
      // saving
      res[t, :m] = aa;
      res[t, (m + 1): ] = symmat_to_vector(PP);
    }
  }
  return res;
}

/**

@section Log-likelihood

*/
/** ssm_lpdf
Log-likelihood of a Linear Gaussian State Space Model

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
@return real The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

The log-likelihood of a linear Gaussian state space model is,
If the the system matrices and initial conditions are known, the log likelihood is
$$
\begin{aligned}[t]
\log L(\mat{Y}_n) &= \log p(\vec{y}_1, \dots, \vec{y}_n) = \sum_{t = 1}^n \log p(\vec{y}_t | \mat{Y}_{t - 1}) \\
&= - \frac{np}{2} \log 2 \pi - \frac{1}{2} \sum_{t = 1}^n \left( \log \left| \mat{F}_t \right| + \vec{v}\T \mat{F}_t^{-1} \vec{v}_t \right)
\end{aligned} ,
$$
where $\mat{F}_t$ and $\mat{V}_t$ come from a forward pass of the Kalman filter.

*/
real ssm_lpdf(vector[] y,
               vector[] d, matrix[] Z, matrix[] H,
               vector[] c, matrix[] T, matrix[] R, matrix[] Q,
               vector a1, matrix P1) {
  real ll;
  int n;
  int m;
  int p;
  int q;
  n = size(y); // number of obs
  m = dims(Z)[2];
  p = dims(Z)[3];
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
    // result matricees for each iteration
    vector[n] ll_obs;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form(Q_t, R_t);

    a = a1;
    P = P1;
    for (t in 1:n) {
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
        if (size(R) > 1 && size(Q) > 1) {
          RQR = quad_form(Q_t, R_t);
        }
      }
      v = ssm_filter_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_filter_update_Finv(P, Z_t, H_t);
      K = ssm_filter_update_K(P, Z_t, T_t, Finv);
      ll_obs[t] = ssm_filter_update_ll(v, Finv);
      // don't save a, P for last iteration
      if (t < n) {
        a = ssm_filter_update_a(a, c_t, T_t, v, K);
        P = ssm_filter_update_P(P, Z_t, T_t, RQR, K);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}

/**

@section Time-Invariant Kalman Filter

*/

/** matrix_diff

@param matrix A An $m \times n$ matrix.
@param matrix B An $m \times n$ matrix.
@return real If converged, then 1, else 0.

The difference between $A$ and $B$ is calculated as,
$$
d(A, B) = \max(A - B) / \max(A)
$$

*/
real matrix_diff(matrix A, matrix B) {
  real eps;
  real norm_AB;
  real norm_A;
  real a;
  real ab;
  int m;
  int n;
  m = rows(A);
  n = cols(A);
  eps = 0.0;
  norm_A = 0.0;
  norm_AB = 0.0;
  for (i in 1:m) {
    for (j in 1:n) {
      a = fabs(A[i, j]);
      ab = fabs(A[i, j] - B[i, j]);
      if (a > norm_A) {
        norm_A = a;
      }
      if (ab > norm_AB) {
        norm_AB = ab;
      }
    }
  }
  eps = norm_AB / norm_A;
  return eps;
}

/** ssm_constant_lpdf
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
@return real The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.

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
    real matdiff;
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
          matdiff = matrix_diff(P, P_old);
          if (matdiff < tol) {
            converged = 1;
          }
        }
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}



/**

@section Common Smoother Functions

*/
/** ssm_smooth_update_r
Update $\vec{r}_t$ in smoothing recursions

@param vector r An $m \times 1$ vector with $\vec{r}_{t-1}$
@param matrix Z A $p \times m$ vector with $\mat{Z}_{t}$
@param vector v A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.
@param matrix Finv A $p \times p$ matrix of the forecast precision, $\mat{F}^{-1}_t$.
@param matrix L An $m \times m$ matrix with $\mat{L}_t$.
@return matrix An $m \times 1$ vector with $\vec{r}_t$.

In smoothing recursions, the vector $\vec{r}_t$ is updated with,
$$
\vec{r}_{t - 1} = \mat{Z}' \mat{F}^{-1}_t \vec{v}_t + \mat{L}' \vec{r}_{t} .
$$

See [@DurbinKoopman2012, Sec 4.4.4, p. 91]
*/
vector ssm_smooth_update_r(vector r, matrix Z, vector v, matrix Finv,
                           matrix L) {
  vector[num_elements(r)] r_new;
  r_new = Z ' * Finv * v + L ' * r;
  return r_new;
}

/** ssm_smooth_update_N
Update $\mat{N}_t$ in smoothing recursions

@param vector N An $m \times 1$ vector with $\vec{N}_{t-1}$
@param matrix Z A $p \times m$ vector with $\mat{Z}_{t}$
@param matrix Finv A $p \times p$ matrix of the forecast precision, $\mat{F}^{-1}_t$.
@param matrix L An $m \times m$ matrix with $\mat{L}_t$.
@return matrix An $m \times m$ matrix with $\vec{N}_t$.

In smoothing recursions, the matrix $\vec{N}_t$ is updated with,
$$
\mat{N}_{t - 1} = \mat{Z}_t' \mat{F}^{-1}_t \mat{Z}_t + \mat{L}_t' \mat{N}_t \mat{L}_t .
$$

See [@DurbinKoopman2012, Sec 4.4.4, p. 91]
*/
matrix ssm_smooth_update_N(matrix N, matrix Z, matrix Finv, matrix L) {
  matrix[rows(N), cols(N)] N_new;
  # may not need this to_symmetric_matrix
  N_new = to_symmetric_matrix(quad_form(Finv, Z) + quad_form(N, L));
  return N_new;
}


/** ssm_smooth_state_size
The number of elements in vectors returned by `ssm_smooth_state`

@param int m The number of states.
@return int The size of the vectors is $m + m (m + 1) / 2$.

*/
int ssm_smooth_state_size(int m) {
  int sz;
  sz = m + symmat_size(m);
  return sz;
}

/** ssm_smooth_state_get_mean
Extract $\hat{\vec{\alpha}}_t$ from vectors returned by `ssm_smooth_state`

@param vector x A vector returned by `ssm_smooth_state`
@param int q The number of state disturbances, $\vec{\eta}_t$.
@return vector An $m \times 1$ vector with $\hat{\vec{\eta}}_t$.
*/
vector ssm_smooth_state_get_mean(vector x, int m) {
  vector[m] alpha;
  alpha = x[ :m];
  return alpha;
}

/** ssm_smooth_state_get_var
Extract $mat{V}_t$ from vectors returned by `ssm_smooth_state`

@param vector x A vector returned by `ssm_smooth_state`
@param int m The number of states
@return matrix An $m \times m$ matrix with $\mat{V}_t$.
*/
matrix ssm_smooth_state_get_var(vector x, int m) {
  matrix[m, m] V;
  V = vector_to_symmat(x[(m + 1): ], m);
  return V;
}


/** ssm_smooth_state
The state smoother

This calculates the mean and variance of the states, $\vec{\alpha}_t$, given the entire sequence, $\vec{y}_{1:n}$.

@param vector[] filter Results of `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@return vector[] An array of vectors constaining $\hat{\vec{\alpha}}_t$ and $\mat{V}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:n})$.
  in the format described below.

For `Z` and `T` the array can have a size of 1, if it is not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `T`) if it is time varying.

The vectors returned by this function have $m + m ^ 2$ elements in this format,
$$
(\hat{\vec{\alpha}}_t', \VEC(\mat{V}_t)' )'.
$$
Use the `ssm_smooth_state_get_mean` and `ssm_smooth_state_get_var` to extract components
from the returned vectors.

value                                        length          start                 end
-------------------------------------------- --------------- --------------------- --------------------
$\hat{\vec{\alpha}}_t$                       $m$             $1$                   $m$
$\mat{V}_t$                                  $m (m + 1) / 2$ $m + 1$               $m + m (m + 1) / 2$


See @DurbinKoopman2012, Eq 4.44 and eq 4.69.
*/
vector[] ssm_smooth_state(vector[] filter, matrix[] Z, matrix[] T) {
  vector[ssm_smooth_state_size(dims(Z)[3])] res[size(filter)];
  int n;
  int m;
  int p;
  n = size(filter);
  m = dims(Z)[3];
  p = dims(Z)[2];
  {
    // system matrices for current iteration
    matrix[p, m] Z_t;
    matrix[m, m] T_t;
    // smoother results
    vector[m] r;
    matrix[m, m] N;
    matrix[m, m] L;
    vector[m] alpha;
    matrix[m, m] V;
    // filter results
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    vector[m] a;
    matrix[m, m] P;

    if (size(Z) == 1) {
      Z_t = Z[1];
    }
    if (size(T) == 1) {
      T_t = T[1];
    }
    // initialize smoother
    // r and N go from n, n - 1, ..., 1, 0.
    // r_n and N_n
    r = rep_vector(0.0, m);
    N = rep_matrix(0.0, m, m);
    // move backwards in time: t, ..., 1
    for (i in 0:(n - 1)) {
      int t;
      t = n - i;
      // set time-varying system matrices
      if (size(Z) > 1) {
        Z_t = Z[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      // get filtered values
      K = ssm_filter_get_K(filter[t], m, p);
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      a = ssm_filter_get_a(filter[t], m, p);
      P = ssm_filter_get_P(filter[t], m, p);
      // updating
      // L_t
      L = ssm_filter_update_L(Z_t, T_t, K);
      // r_{t - 1} and N_{t - 1}
      r = ssm_smooth_update_r(r, Z_t, v, Finv, L);
      N = ssm_smooth_update_N(N, Z_t, Finv, L);
      // hat(alpha)_{t} and V_t which use r and N from (t - 1)
      alpha = a + P * r;
      V = to_symmetric_matrix(P - P * N * P);
      // saving
      res[t, :m] = alpha;
      res[t, (m + 1): ] = symmat_to_vector(V);
    }
  }
  return res;
}


/** ssm_smooth_eps_size
The size of the vectors returned by `ssm_smooth_eps`

@param int p The length of the observation vectors, $\vec{y}_t$.
@return int The size of the vectors is $p + p (p + 1) / 2$.

*/
int ssm_smooth_eps_size(int p) {
  int sz;
  sz = p + symmat_size(p);
  return sz;
}

/** ssm_smooth_eps_get_mean

Extract $\hat{\vec{\varepsilon}}_t$ from vectors returned by `ssm_smooth_eps`

@param x A vector from the results of `ssm_smooth_eps`.
@param int p The length of the observation vectors, $\vec{y}_t$.
@return vector A $p \times 1$ vector with $\hat{\vec{\varepsilon}}_t$.
*/
vector ssm_smooth_eps_get_mean(vector x, int p) {
  vector[p] eps;
  eps = x[ :p];
  return eps;
}

/** ssm_smooth_eps_get_var
Extract $\Var(\varepsilon_t|\vec{y}_{1:n})$ from vectors returned by `ssm_smooth_eps`

@param vector x A vector returned by `ssm_smooth_eps`
@param int p The length of the observation vectors, $\vec{y}_t$.
@return matrix A $p \times p$ matrix with $\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$

*/
matrix ssm_smooth_eps_get_var(vector x, int p) {
  matrix[p, p] eps_var;
  eps_var = vector_to_symmat(x[(p + 1): ], p);
  return eps_var;
}

/** ssm_smooth_eps

The observation disturbance smoother

This calculates the mean and variance of the observation disturbances, $\vec{\varepsilon}_t$,
given the entire sequence, $\vec{y}_{1:n}$.


@param vector[] filter Results of `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@return vector[] An array of vectors constaining $\hat{\vec{\varepsilon}}_t$ and $\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$ in the format described below.


For Z`, `H`, T`, the array can have a size of 1, if it is not time-varying, or a size of $n$ (for `Z`, `H`) or $n - 1$ (for `T`),
if it is time varying.

The vectors returned by this function have $p + p (p + 1) / 2$ elements in this format,
$$
(\hat{\vec{\varepsilon}}_t', \VEC(\Var(\vec{\varepsilon}_t | \vec{y}_{1:n}))' )'
$$

value                                        length          start                 end
-------------------------------------------- --------------- --------------------- --------------------
$\hat{\vec{\varepsilon}}_t$                  $p$             $1$                   $p$
$\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$  $p (p + 1) / 2$ $p + 1$               $p + p (p + 1) / 2$


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
    // r and N go from n, n - 1, ..., 1, 0.
    // r_n and N_n
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
      // r_{t - 1} and N_{t - 1}
      r = ssm_smooth_update_r(r, Z_t, v, Finv, L);
      N = ssm_smooth_update_N(N, Z_t, Finv, L);
      // eps_t and V(eps_t|y)
      eps = H_t * (Finv * v - K ' * r);
      var_eps = to_symmetric_matrix(H_t - H_t * (Finv + quad_form(N, K)) * H_t);
      // saving
      res[t, :p] = eps;
      res[t, (p + 1): ] = symmat_to_vector(var_eps);
    }
  }
  return res;
}

/** ssm_smooth_eta
The size of the vectors returned by `ssm_smooth_eta`

@param int p The length of the observation vectors, $\vec{y}_t$.
@return int The size of the vectors is $q + q (q + 1) / 2$.

*/
int ssm_smooth_eta_size(int q) {
  int sz;
  sz = q + symmat_size(q);
  return sz;
}

/** ssm_smooth_eta_get_mean
Extract $\hat{\vec{\varepsilon}}_t$ from vectors returned by `ssm_smooth_eta`

@param vector x A vector returned by `ssm_smooth_eta`
@param int q The number of state disturbances, $\vec{\eta}_t$.
@return vector A $q \times 1$ vector with $\hat{\vec{\eta}}_t$.

*/
vector ssm_smooth_eta_get_mean(vector x, int q) {
  vector[q] eta;
  eta = x[ :q];
  return eta;
}

/** ssm_smooth_eta_get_var
Extract $\Var(\eta_t|\vec{y}_{1:n})$ from vectors returned by `ssm_smooth_eta`

@param vector x A vector returned by `ssm_smooth_eta`
@param int q The number of state disturbances, $\vec{\eta}_t$.
@return matrix A $q \times q$ matrix with $\Var(\vec{\eta}_t | \vec{y}_{1:n})$.

*/
matrix ssm_smooth_eta_get_var(vector x, int q) {
  matrix[q, q] eta_var;
  eta_var = vector_to_symmat(x[(q + 1): ], q);
  return eta_var;
}

/** ssm_smooth_eta

The state disturbance smoother

This calculates the mean and variance of the observation disturbances, $\vec{\eta}_t$,
given the entire sequence, $\vec{y}_{1:n}$.

@param vector[] filter Results of `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@return vector[] An array of vectors constaining $\hat{\vec{\eta}}_t$ and $\Var(\vec{\eta}_t | \vec{y}_{1:n})$ in the format described below.

For `Z`, `T`, `R`, `Q` the array can have a size of 1, if it is not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `T`, `R`, `Q`) if it is time varying.

The vectors returned by this function have $q + q (q + 1) / 2$ elements in this format,
$$
(\hat{\vec{\eta}}_t', \VEC(\Var(\vec{\eta}_t | \vec{y}_{1:n}))' ).
$$
Use the `ssm_smooth_eta_get_mean` and `ssm_smooth_eta_get_var` to extract components
from the returned vectors.

value                                 length          start                 end
------------------------------------- --------------- --------------------- --------------------
$\hat{\vec{\eta}}_t$                  $q$             $1$                   $q$
$\Var(\vec{\eta}_t | \vec{y}_{1:n})$  $q (q + 1) / 2$ $q + 1$               $q + q (q + 1) / 2$

See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]

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
    for (i in 0:(n - 1)) {
      int t;
      // move backwards in time
      t = n - i;
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
      res[t, :q] = eta;
      res[t, (q + 1): ] = symmat_to_vector(var_eta);
    }
  }
  return res;
}
/** ssm_smooth_faststate

The fast state smoother

The fast state smoother calculates $\hat{\vec{\alpha}}_t = \E(\vec{\alpha}_t | \vec{y}_{1:n})$.
$$
\hat{\vec{\alpha}}_{t + 1} = \mat{T}_t \hat{\vec{\alpha}}_{t} + \mat{R}_t \mat{Q}_t \mat{R}'_t \vec{r}_t ,
$$
where $r_t$ is calcualted from the state disturbance smoother.
The smoother is initialized at $t = 1$ with $\hat{\vec{\alpha}}_t = \vec{a}_1 + \mat{P}_1 \vec{r}_0$.

Unlike the normal state smoother, it does not calculate the variances of the smoothed state.

@param vector[] filter The results of `ssm_filter`
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param vector[] c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@return vector[] An array of size $n$ of $m \times 1$ vectors containing $\hat{\vec{\alpha}}_t$.

For  `Z`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]

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
    vector[m] r[n + 1];
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
    // r goes from t = n, ..., 1, 0.
    // r_n
    r[n + 1] = rep_vector(0.0, m);
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
      // r_{t - 1}
      r[t] = ssm_smooth_update_r(r[t + 1], Z_t, v, Finv, L);
    }
    // calculate smoothed states
    a1 = ssm_filter_get_a(filter[1], m, p);
    P1 = ssm_filter_get_P(filter[1], m, p);
    // r[1] = r_0
    alpha[1] = a1 + P1 * r[1];
    // 1:(n - 1) -> \alpha_{2}:\alpha_{n}
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
      // `r[t + 1]` = $r_{t}$
      // alpha_{t + 1} = c_t + T_t * \alpha_t + R_t Q_t R'_t r_t
      alpha[t + 1] = c_t + T_t * alpha[t] + RQR * r[t + 1];
    }
  }
  return alpha;
}




/**
@section Simulators and Smoothing Simulators

*/
/** ssm_sim_idx

Indexes of each component of `ssm_sim_rng` results.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return int[,] A 4 x 3 array of integers

The returned array has columns (length, start location, and end location)
for rows: $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\varepsilon}_t$, and $\vec{\eta}_t$ in the results of `ssm_sim_rng`.

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

/** ssm_sim_size

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

/** ssm_sim_get_y

Extract $\vec{y}_t$ from vectors returned by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector vector A $p \times 1$ vector with $\vec{y}_t$.

*/
vector ssm_sim_get_y(vector x, int m, int p, int q) {
  vector[p] y;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  y = x[idx[1, 2]:idx[1, 3]];
  return y;
}

/** ssm_sim_get_a

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

/** ssm_sim_get_eps

Extract $\vec{\varepsilon}_t$ from vectors returne by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector vector A $p \times 1$ vector with $\vec{\varepsilon}_t$.


*/
vector ssm_sim_get_eps(vector x, int m, int p, int q) {
  vector[p] eps;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  eps = x[idx[3, 2]:idx[3, 3]];
  return eps;
}

/** ssm_sim_get_eta
Extract $\vec{\eta}_t$ from vectors returne by `ssm_sim_rng`.

@param int m The number of states
@param int p The length of the observation vector
@param int q The number of state disturbances
@return vector vector A $q \times 1$ vector with $\vec{\eta}_t$.

*/
vector ssm_sim_get_eta(vector x, int m, int p, int q) {
  vector[q] eta;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  eta = x[idx[4, 2]:idx[4, 3]];
  return eta;
}

/** ssm_sim_rng

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
@return Array of size $n$ of vectors with Draw $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\eta}_t$ and $\vec{\varepsilon}_t$. See the description.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

Draw $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\eta}_t$ and $\vec{\varepsilon}_t$ from
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

element         length         start         end
--------------- -------------- ------------- -----------
$y_t$           $p$            $1$           $p$
$\alpha$_t      $m$            $p + 1$       $p + m$
$\varepsilon_t$ $p$            $p + m + 1$   $2 p + m$
$\eta_t$        $q$            $2 p + m + 1$ $2 p + m + q$

It is preferrable to use `ssm_sim_get_y`, `ssm_sim_get_a`, `ssm_sim_get_eps`,
and `ssm_sim_get_eta` to extract values from these vectors.

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

/**

@section Simulation Smoothers

*/

/** ssm_simsmo_state_rng
State simulation smoother

Draw samples from the posterior distribution of the states,
$\tilde{\vec{\alpha}}_{1:n} \sim p(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.

@param vector[] alpha An of size $n$ of $m \times 1$ vectors containing the smoothed expected values of the states, $\E(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$. These are returned by `sim_smooth_faststates`. If `sim_smooth_state` was used, then the expected values need to first be extracted using `sim_smooth_state_get_mean`.
@param vector[] d Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param vector[] c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@param vector a1 Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
@param matrix P1 Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
@return vector[] Array of size $n$ of $m \times 1$ vectors containing a single draw from $(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].

*/
vector[] ssm_simsmo_states_rng(vector[] alpha,
                      vector[] d, matrix[] Z, matrix[] H,
                      vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                      vector a1, matrix P1) {
    vector[dims(Z)[2]] draws[size(alpha)];
    int n;
    int p;
    int m;
    int q;
    n = size(alpha);
    p = dims(Z)[2];
    m = dims(Z)[3];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[p] y[n];
      vector[m] alpha_hat_plus[n];
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter with simulated y's
      filter = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      alpha_hat_plus = ssm_smooth_faststate(filter, c, Z, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha[i]);
      }
    }
    return draws;
}

/** ssm_simsmo_eta_rng

State disturbance simulation smoother

Draw samples from the posterior distribution of the observation disturbances,
$\tilde{\vec{\eta}}_{1:n} \sim p(\vec{\eta}_{1:n} | \vec{y}_{1:n})$.

@param vector[] eta Values returned by `sim_smooth_eta`
@param vector[] d Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param vector[] c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@param vector a1 Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
@param matrix P1 Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
@return vector[] Array of size $n$ of $q \times 1$ vectors containing a single draw from $(\vec{\eta}_{1:n} | \vec{y}_{1:n})$.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].

*/
vector[] ssm_simsmo_eta_rng(vector[] eta,
                            vector[] d, matrix[] Z, matrix[] H,
                            vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                            vector a1, matrix P1) {
    vector[dims(Q)[2]] draws[size(eta)];
    int n;
    int p;
    int m;
    int q;
    n = size(eta);
    p = dims(Z)[2];
    m = dims(Z)[3];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter[n];
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eta_size(q)] etahat_plus[n];
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct eta samples
      etahat_plus = ssm_smooth_eta(filter, Z, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_eta(sims[i], m, p, q)
                                    - ssm_smooth_eta_get_mean(etahat_plus[i], q)
                                    + ssm_smooth_eta_get_mean(eta[i], q));
      }
    }
    return draws;
}

/** ssm_simsmo_eps_rng

Observation disturbance simulation smoother

Draw samples from the posterior distribution of the observation disturbances,
$\tilde{\vec{\varepsilon}}_{1:n} \sim p(\vec{\varepsilon}_{1:n} | \vec{y}_{1:n})$.

@param vector[] eps Values returned by `sim_smooth_eps`
@param vector[] d Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
@param matrix[] Z Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
@param matrix[] H Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
@param vector[] c State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
@param matrix[] T Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
@param matrix[] R State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
@param matrix[] Q State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
@param vector a1 Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
@param matrix P1 Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
@return vector[] Array of size $n$ of $p \times 1$ vectors containing a single draw from $(\vec{\varepsilon}_{1:n} | \vec{y}_{1:n})$.

For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].

*/
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

/**

@section Stationary

*/
/** pacf_to_acf

Partial Autocorrelations to Autocorrelations

@param vector x A vector of coefficients of a partial autocorrelation function
@return vector A vector of coefficients of an Autocorrelation function

*/
vector pacf_to_acf(vector x) {
  matrix[num_elements(x), num_elements(x)] y;
  int n;
  n = num_elements(x);
  y = rep_matrix(0.0, n, n);
  for (k in 1:n) {
    for (i in 1:(k - 1)) {
      y[k, i] = y[k - 1, i] + x[k] * y[k - 1, k - i];
    }
    y[k, k] = x[k];
    print(y);
  }
  return -y[n] ';
}

/** constrain_stationary
Constrain vector of coefficients to the stationary and intertible region for AR or MA functions.

@param vector x An unconstrained vector in $(-\infty, \infty)$
@return vector A vector of coefficients for a stationary AR or inverible MA process.

See @Jones1980a, @Jones1987a, @Monahan1984a, @AnsleyKohn1986a, and the functions
`tools.constrain_stationary_univariate` and `tools.unconstraine_stationary_univariate` in
[statsmodels.tsa.statespace](http://www.statsmodels.org/dev/statespace.html#statespace-tools).

*/
vector constrain_stationary(vector x) {
  vector[num_elements(x)] r;
  int n;
  n = num_elements(x);
  // transform (-Inf, Inf) to (-1, 1)
  for (i in 1:n) {
    r[i] = x[i] / (sqrt(1.0 + pow(x[i], 2)));
  }
  // Transform PACF to ACF
  return pacf_to_acf(r);
}

/** acf_to_pacf
Convert coefficients of an autocorrelation function to partial autocorrelations.

@param vector x Coeffcients of an autocorrelation function.
@return vector A vector of coefficients of the corresponding partial autocorrelation function.

*/
vector acf_to_pacf(vector x) {
  matrix[num_elements(x), num_elements(x)] y;
  vector[num_elements(x)] r;
  int n;
  n = num_elements(x);
  y = rep_matrix(0.0, n, n);
  y[n] = -x ';
  for (j in 0:(n - 1)) {
    int k;
    k = n - j;
    for (i in 1:(k - 1)) {
      y[k - 1, i] = (y[k, i] - y[k, k] * y[k, k - i]) / (1 - pow(y[k, k], 2));
    }
  }
  r = diagonal(y);
  return r;
}

/** unconstrain_stationary
Transform from stationary and invertible space to $(-\infty, \infty)$.

@param vector x Coeffcients of an autocorrelation function.
@return vector Coefficients of the corresponding partial autocorrelation function.

*/
vector unconstrain_stationary(vector x) {
  matrix[num_elements(x), num_elements(x)] y;
  vector[num_elements(x)] r;
  vector[num_elements(x)] z;
  int n;
  n = num_elements(x);
  // Transform ACF to PACF
  r = acf_to_pacf(x);
  // Transform (-1, 1) to (-Inf, Inf)
  for (i in 1:n) {
    z[i] = r[i] / (sqrt(1.0 - pow(r[i], 2)));
  }
  return z;
}

/** kronecker_prod

Kronecker product

The Kronecker product of a $A$ and $B$ is
$$
A \otimes B =
\begin{bmatrix}
a_{11} B \cdots a_{1n} B \\
\vdots & \ddots & vdots \\
a_{m1} B & \cdots & a_{mn} B
\end{bmatrix} .
$$

@param matrix A An $m \times n$ matrix
@param matrix B A $p \times q$ matrix
@return matrix An $mp \times nq$ matrix.

*/
matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}

/** arima_stationary_cov

Find the covariance of the stationary distribution of an ARMA model

@param matrix T The $m \times m$ transition matrix
@param matrix R The $m \times q$ system disturbance selection matrix
@return matrix An $m \times m$ matrix with the stationary covariance matrix.

The initial conditions are $\alpha_1 \sim N(0, \sigma^2 Q_0)$,
where $Q_0$ is the solution to
$$
(T \otimes T) \VEC(Q_0) = \VEC(R R')
$$
where $\VEC(Q_0)$ and $\VEC(R R')$ are the stacked columns of $Q_0$ and $R R'$

See @DurbinKoopman2012, Sec 5.6.2.

*/
matrix arima_stationary_cov(matrix T, matrix R) {
  matrix[rows(T), cols(T)] Q0;
  matrix[rows(T) * rows(T), rows(T) * rows(T)] TT;
  vector[rows(T) * rows(T)] RR;
  int m;
  int m2;
  m = rows(T);
  m2 = m * m;
  RR = to_vector(tcrossprod(R));
  TT = kronecker_prod(T, T);
  Q0 = to_matrix_colwise((diag_matrix(rep_vector(1.0, m2)) - TT) \ RR, m, m);
  return Q0;
}
