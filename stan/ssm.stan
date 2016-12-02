/**
`r render_chapter("Stan Functions")`

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

`r render_section("Utility Functions")`

*/
/**
---
function: to_symmetric_matrix
args:
  - name: x
    description: An $n \times n$ matrix
returns: An $n \times n$ symmetric matrix, $0.5 (x + x')$.
---
Ensure a matrix is symmetric
*/
matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}

/**
---
function: to_matrix_colwise
args:
  - name: v
    description: An $n \times m$ vector.
  - name: m
    description: Number of rows in the vector
  - name: n
    description: Number of columns in the vector
returns: A $m \times n$ matrix containting the elements from `v`
---
Convert vector to a matrix (column-major).
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

/**
---
function: matrix_pow
args:
- name: A
  description: The matrix to take the power of
- name: n
  description: The order of the power. This is sepcified as a real number to avoid compiler warnings, but only the integer part is used.
---
Calculate the power of a matrix, $\mat{A}^n$.

*/
matrix matrix_pow(matrix A, real n);


matrix matrix_pow(matrix A, real n) {
  real nn;
  nn = floor(n);
  if (nn == 0) {
    return diag_matrix(rep_vector(1., rows(A)));
  } else if (nn == 1) {
    return A;
  } else if (nn > 1) {
    # recurively this is n log n.
    if (fmod(nn, 2.) > 0) {
      # If odd
      return A * matrix_pow(A, nn - 1);
    } else {
      # If even
      return matrix_pow(A, nn / 2) * matrix_pow(A, nn / 2);
    }
  } else {
    # n < 0
    reject("Only non-negative values of n are allowed");
    return A;
  }
}

// // This is the definition in the Stan book and uses an int n.
// // It requires n multiplications
//
// matrix matrix_pow(matrix A, int n);
//
// matrix matrix_pow(matrix A, int n) {
//   if (n == 0.) {
//     return diag_matrix(rep_vector(1., rows(A)));
//   } else {
//     return a * matrix_pow(A, n - 1);
//   }


/**
---
function: symmat_size
args:
- name: n
  description: The number of rows and columns in the matrix.
returns: The number of unique elements
---

Calculate the number of unique elements in a symmetric matrix

The number of unique elements in an $m \times m$ matrix is
$(m \times (m + 1)) / 2$.


*/

int symmat_size(int n) {
  int sz;
  // This is calculated iteratively to avoid the Stan warning for
  // integer division
  sz = 0;
  for (i in 1:n) {
    sz = sz + i;
  }
  return sz;
}

/**
---
function: find_symmat_dim
args:
- name: n
  description: The number of unique elements in a symmetric matrix.
returns: The dimension of the associated symmetric matrix.
---

Given vector with $n$ elements containing the $m (m + 1) / 2$ elements of a symmetric matrix,
return $m$.


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

/**
---
function: vector_to_symmat
args:
- name: x
  description: The vector with the unique elements
- name: n
  description: The dimensions of the returned matrix, $n \times n$.
returns: matrix An $n \times n$ symmetric matrix.
---
Convert a vector to a symmetric matrix

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

/**
---
function: symmat_to_vector
args:
- name: x
  description: An $n \times n$ matrix.
returns: A $n (n + 1) / 2$ vector with the unique elements in $x$.
---


Convert an $n \times n$ symmetric matrix to a length $n (n + 1) / 2$ vector
containing its unique elements.


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
---
function: rep_lower_triangular_matrix
args:
- name: x
  description: Value used for the non-zero elements of the matrix.
- name: m
  description: number of rows
- name: n
  description: number of columns
- name: diag
  description: If true, then include 1's on the diagonal.
returns: An $m \times n$ lower triangular matrix
---


Fill in an lower triangular matrix.

*/

matrix rep_lower_triangular_matrix(real x, int m, int n, int diag) {
  matrix[m, n] A;
  for (i in 1:m) {
    for (j in 1:n) {
      if (i > j) {
        A[i, j] = x;
      } else if (i == j) {
        if (diag) {
          A[i, j] = x;
        } else {
          A[i, j] = 0.;
        }
      } else {
        A[i, j] = 0.;
      }
    }
  }
  return A;
}

/**
---
function: rep_upper_triangular_matrix
args:
- name: x
  description: Value used for the non-zero elements of the matrix.
- name: m
  description: number of rows
- name: n
  description: number of columns
- name: diag
  description: If true, then include 1's on the diagonal.
returns: An $m \times n$ upper triangular matrix
---


Fill in an upper triangular matrix

*/

matrix rep_upper_triangular_matrix(real x, int m, int n, int diag) {
  matrix[m, n] A;
  for (i in 1:m) {
    for (j in 1:n) {
      # if row less than column
      if (i < j) {
        A[i, j] = x;
      } else if (i == j) {
        if (diag) {
          A[i, j] = x;
        } else {
          A[i, j] = 0.;
        }
      } else {
        A[i, j] = 0.;
      }
    }
  }
  return A;
}


/**
---
function: rep_diagonal_matrix
args:
- name: x
  description: Value used for the non-zero elements of the matrix.
- name: m
  description: number of rows
- name: n
  description: number of columns
- name: k
  description: Index of the diagonal
returns: An $m \times n$ upper triangular matrix
---


Create a diagonal $m \times n$ matrix with values $x$ on the $k$-th diagonal.

*/

matrix rep_diagonal_matrix(real x, int m, int n, int k) {
  matrix[m, n] A;
  int mn;
  A = rep_matrix(0., m, n);
  mn = min(m, n);
  if (k >= 0) {
    for (i in 1:min(m, n - k)) {
      A[i, i + k] = x;
    }
  } else {
    for (i in 1:min(m + k, n)) {
      A[i - k, i] = x;
    }
  }
  return A;
}

/**
---
function: fill_matrix
args:
- name: x
  description: A $p \times q$, $p \leq m$, $\q \leq n$ matrix
- name: m
  description: Number of rows in the returned matrix
- name: n
  description: Number of columns in the returned matrix
- name: i
  description: Indices mapping the rows of $A$ to the rows in the output matrix
- name: j
  description: Indices mapping the columns of $A$ to the columns of the output matrix
- name: a
  description: The default value in the returned matrix
returns: A $m \times n$ matrix
---


Given a $p \times q$ matrix $\mat{X}$, default value $a$, and indexes $\vec{I} = i_1, ..., i_p$,
and $\vec{J} = j_1, ...j_q$, return a $m \times n$ matrix where $m \geq p$, $n \geq q$, where
$$
Y_{k, l} =
\begin{cases}
X_{i, j} & \text{if $k = i$, $l = j$, for some $i \in \vec{I}$, $j \in \vec{J}$,} \\
a & \text{otherwise} .
\end{cases}
$$


*/

matrix fill_matrix(matrix x, int m, int n, int[] i, int[] j, real a) {
  matrix[m, n] ret;
  ret = rep_matrix(a, m, n);
  ret[i, j] = x;
  return ret;
}

/**
---
function: fill_vector
args:
- name: x
  description: A $p \times q$, $p \leq m$, $\q \leq n$ matrix
- name: n
  description: Number of elements in the returned vector
- name: i
  description: Indices mapping the rows of $A$ to the rows in the output matrix
- name: a
  description: The default value in the returned vector
returns: A $n \times 1$ matrix
---


Given an $m \times 1$ vector $\vec{x}$, an integer $n \geq m$, a default value $a$,
and indexes $\vec{I} = i_1, ..., i_m \in 1:n$, return a $n \times 1$ vector where
y_{j} =
\begin{cases}
x_{i} & \text{if $j = i$ for some $i \in \vec{I}$,} \\
a & \text{otherwise}
\end{cases} .
$$

*/

vector fill_vector(vector x, int n, int[] i, real a) {
  vector[n] ret;
  ret = rep_vector(a, n);
  ret[i] = x;
  return ret;
}

/**
---
function: int_sum_true
args:
- name: x
  description: An array of length $n$ of integers
returns: An integer between 0 and $n$.
---


For an array of integers, return the indexes where it is greater than zero.


*/

int int_sum_true(int[] x) {
  int n;
  n = 0;
  for (i in 1:num_elements(x)) {
    if (int_step(x[i])) {
      n = n + 1;
    }
  }
  return n;
}

/**
---
function: int_sum_false
args:
- name: x
  description: An array of length $n$ of integers
returns: An integer between 0 and $n$.
---


For an array of integers, return the indexes where it is less than or equal to zero.


*/

int int_sum_false(int[] x) {
  int n;
  n = 0;
  for (i in 1:num_elements(x)) {
    if (! int_step(x[i])) {
      n = n + 1;
    }
  }
  return n;
}


/**
---
function: mask_indexes
args:
- name: x
  description: An array of length $n$ of integers
- name: n
  description: The number of false values in `x`.
returns: An array of integers with elements having values between 1 and $m$.
---


For an array of integers, `x`, return the indexes where
mask is not true (`x[i] <= 0`).
The primary use of this function is where `x` represents
indicators for missing values,  and it is used to extract
the indexes of non-missing values.


*/

int[] mask_indexes(int[] x, int n) {
  int idx[n];
  int j;
  j = 1;
  if (n > 0) {
    for (i in 1:num_elements(x)) {
      if (! int_step(x[i]) && j <= n) {
        idx[j] = i;
        j = j + 1;
      }
    }
  }
  return idx;
}


/**
---
function: select_indexes
args:
- name: x
  description: An array of length $m$ of integers
- name: n
  description: The number of true values in `x`.
returns: An array of integers with elements having values between 1 and $m$.
---


For an array of integers, `x`, return the indexes where
the elements are true (`x[i] > 0`).
The primary use of this function is where `x` represents
indicators for non-missing values, and it is used to extract
the indexes of non-missing values.


*/

int[] select_indexes(int[] x, int n) {
  int idx[n];
  int j;
  j = 1;
  if (n > 0) {
    for (i in 1:num_elements(x)) {
      if (int_step(x[i]) && j <= n) {
        idx[j] = i;
        j = j + 1;
      }
    }
  }
  return idx;
}

/**
---
function: normal2_rng
args:
- name: mu
  description: mean
- name: sigma
  description: variance
returns: A value drawn from the specified normal distribution.
---


Draw samples from a normal distribution with mean `mu` and scale `sigma`.
Unlike the built-in `normal_rng()`, this allows for `sigma = 0`.


*/

real normal2_rng(real mu, real sigma) {
  real y;
  if (sigma <= 0) {
    y = mu;
  } else {
    y = normal_rng(mu, sigma);
  }
  return y;
}

/**
---
function: cholesky_decompose2
args:
- name: A
  description: An $n \times n$ matrix
returns: An $n \times n$ lower-triangular matrix
---


Calculate the Cholesky decomposition of a matrix. Unlike the built-in
function, this handles cases in which the matrix has 0's on the diagonal.

*/

matrix cholesky_decompose2(matrix A) {
  matrix[rows(A), cols(A)] L;
  int n;
  int nonzero[rows(A)];
  int num_nonzero;
  n = rows(A);
  for (i in 1:n) {
    nonzero[i] = (A[i, i] > 0);
  }
  num_nonzero = sum(nonzero);
  if (num_nonzero == n) {
    L = cholesky_decompose(A);
  } else if (num_nonzero == 0) {
    L = rep_matrix(0.0, n, n);
  } else {
    int idx[num_nonzero];
    vector[n] eps;
    idx = select_indexes(nonzero, num_nonzero);
    L = rep_matrix(0.0, n, n);
    L[idx, idx] = cholesky_decompose(A[idx, idx]);
  }
  return L;
}


/**
---
function: multi_normal2_rng
args:
- name: mu
  description: An $n \times 1$ vector of the means
- name: Sigma
  description: An $n \times n$ lower triangular matrix with covariance matrix.
returns: An $n \times 1$ vector drawn from the specified multivariate normal distribution.
---


Sample from a multivariate normal distribution.
Unlike the built-in `multi_normal_rng`,
this function will still draw samples for deterministic elements in the vector.


*/

vector multi_normal2_rng(vector mu, matrix Sigma) {
  vector[num_elements(mu)] y;
  int n;
  int nonzero[num_elements(mu)];
  int num_nonzero;
  n = num_elements(mu);
  for (i in 1:n) {
    nonzero[i] = (Sigma[i, i] > 0);
  }
  num_nonzero = sum(nonzero);
  if (num_nonzero == n) {
    y = multi_normal_rng(mu, Sigma);
  } else if (num_nonzero == 0) {
    y = mu;
  } else {
    int idx[num_nonzero];
    vector[n] eps;
    idx = select_indexes(nonzero, num_nonzero);
    eps = rep_vector(0.0, n);
    eps[idx] = multi_normal_rng(rep_vector(0.0, num_nonzero), Sigma[idx, idx]);
    y = mu + eps;
  }
  return y;
}

/**
---
function: multi_normal_cholesky2_rng
args:
- name: mu
  description: An $n \times 1$ vector of the means
- name: L
  description: An $n \times n$ lower triangular matrix with the Cholesky decomposition of the covariance matrix.
returns: An $n \times 1$ vector drawn from the specified multivariate normal distribution.
---


Sample from a multivariate normal distribution, parameterized with the Cholesky
decomposition of the covariance matrix. Unlike the built-in `multi_normal_cholesky_rng`,
this function will still draw samples for deterministic elements in the vector.


*/

vector multi_normal_cholesky2_rng(vector mu, matrix L) {
  vector[num_elements(mu)] y;
  int n;
  int nonzero[num_elements(mu)];
  int num_nonzero;
  n = num_elements(mu);
  for (i in 1:n) {
    nonzero[i] = (L[i, i] > 0);
  }
  num_nonzero = sum(nonzero);
  if (num_nonzero == n) {
    y = multi_normal_cholesky_rng(mu, L);
  } else if (num_nonzero == 0) {
    y = mu;
  } else {
    int idx[num_nonzero];
    vector[n] eps;
    idx = select_indexes(nonzero, num_nonzero);
    eps = rep_vector(0.0, n);
    eps[idx] = multi_normal_cholesky_rng(rep_vector(0.0, num_nonzero),
                                         L[idx, idx]);
    y = mu + eps;
  }
  return y;
}


/**

`r render_section("Filtering")`

Functions used in filtering and log-likelihood calculations.
*/
/**
---
function: ssm_update_a
args:
- name: a
  description: An $m \times 1$ vector with the predicted state, $\vec{a}_t$.
- name: c
  description: An $m \times 1$ vector with the system intercept, $\vec{c}_t$
- name: T
  description: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- name: v
  description: A $p \times 1$ vector with the forecast error, $\vec{v}_t$.
- name: K
  description: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
returns: A $m \times 1$ vector with the predicted state at $t + 1$, $\vec{a}_{t + 1}$.
---

Update the expected value of the predicted state, $\vec{a}_{t + 1} = \E(\vec{\alpha}_{t + 1} | \vec{y}_{1:t})$,


The predicted state $\vec{a}_{t + 1}$ is,
$$
\vec{a}_{t + 1} = \mat{T}_t \vec{a}_t + \mat{K}_t \vec{v}_t + \vec{c}_t .
$$

*/

vector ssm_update_a(vector a, vector c, matrix T, vector v, matrix K) {
  vector[num_elements(a)] a_new;
  a_new = T * a + K * v + c;
  return a_new;
}

/**
---
function: ssm_update_P
args:
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: T
  description: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- name: RQR
  description: A $m \times m$ matrix with the system covariance matrix, $\mat{R}_t \mat{Q}_t \mat{R}_t'$.
- name: K
  description: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
returns: An $m \times m$ matrix with the variance of the state, $\vec{P}_{t + 1}$.
---

Update the variance of the state in $t + 1$, $\mat{P}_{t + 1} = \Var(\alpha_{t + 1} | \vec{y}_{1:t})$,

The predicted state variance $\mat{P}_{t + 1}$ is,
$$
\mat{P}_{t + 1} = \mat{T}_t \mat{P}_t (\mat{T}_t - \mat{K}_t \mat{Z}_t)' + \mat{R}_t \mat{Q}_t \mat{R}_t' .
$$

*/

matrix ssm_update_P(matrix P, matrix Z, matrix T,
                           matrix RQR, matrix K) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(T * P * (T - K * Z)' + RQR);
  return P_new;
}

/**
---
function: ssm_update_v
args:
- name: y
  description: A $p \times 1$ vector of the observations, $\vec{y}_t$.
- name: a
  description: A $m \times 1$ vector of the states, $\vec{a}_t$.
- name: d
  description: A $p \times 1$ vector with the observation intercept, $\vec{d}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
returns: A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.
---

Update the forcast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}})$


The forecast error $\vec{v}_t$ is
$$
\vec{v}_t =\vec{y}_t - \mat{Z}_t \vec{a}_t - \vec{d}_t .
$$

*/

vector ssm_update_v(vector y, vector a, vector d, matrix Z) {
  vector[num_elements(y)] v;
  v = y - Z * a - d;
  return v;
}

/**
---
function: ssm_update_v_miss
args:
- name: y
  description: A $p \times 1$ vector of the observations, $\vec{y}_t$.
- name: a
  description: A $m \times 1$ vector of the states, $\vec{a}_t$.
- name: d
  description: A $p \times 1$ vector with the observation intercept, $\vec{d}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: p_t
  description: The number of non-missing elements of the observation vector $\vec{y}_t$.
- name: y_idx
  description: A length $p$ array of integers indexes with the indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.
returns: A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.
---


Update the forcast error, but unlike `ssm_update_v`, allow for missing
values in $\vec{y}_t$.


The elements of the forecast error $\vec{v}_t$ is
$$
\vec{v}_t =
\begin{cases}
  y_{j,t} - \vec{Z}_{j,.,t} \vec{a}_t - d_{j,t} & \text{if $y_{j,t} not missing.} \\
  0 & \text{if $y_{j,t}$ is missing.}
\end{cases}
$$

*/

vector ssm_update_v_miss(vector y, vector a, vector d, matrix Z,
                                int p_t, int[] y_idx) {
  vector[num_elements(y)] v;
  int p;
  p = num_elements(y);
  if (p_t < p) {
    v = rep_vector(0., p);
    if (p_t > 0) {
      int idx[p_t];
      vector[p_t] y_star;
      vector[p_t] d_star;
      matrix[p_t, cols(Z)] Z_star;
      idx = y_idx[1:p_t];
      y_star = y[idx];
      d_star = d[idx];
      Z_star = Z[idx, :];
      v[idx] = ssm_update_v(y_star, a, d_star, Z_star);
    }
  } else {
    v = ssm_update_v(y, a, d, Z);
  }
  return v;
}

/**
---
function: ssm_update_F
args:
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: H
  description: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.
returns: A $p \times p$ vector with $\mat{F}_t$.
---

Update the variance of the forcast error, $\mat{F}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))$


The variance of the forecast error $\mat{F}_t$ is
$$
\mat{F}_t = \mat{Z}_t \mat{P}_t \mat{Z}_t + \mat{H}_t .
$$

*/

matrix ssm_update_F(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] F;
  F = to_symmetric_matrix(quad_form(P, Z') + H);
  return F;
}

/**
---
function: ssm_update_Finv
args:
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: H
  description: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.
returns: A $p \times p$ vector with $\mat{F}^{-1}_t$.
---

Update the precision of the forcast error, $\mat{F}^{-1}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))^{-1}$


This is the inverse of $\mat{F}_t$.

*/

matrix ssm_update_Finv(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] Finv;
  // if can guarantee that F is spd, then take spd inverse.
  Finv = inverse_spd(to_symmetric_matrix(quad_form(P, Z') + H));
  // Finv = inverse(quad_form(P, Z') + H);
  return Finv;
}

/**
---
function: ssm_update_Finv_miss
args:
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: H
  description: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.
- name: p_t
  description: The number of non-missing elements in the observation vector, $\vec{y}_t$.
- name: y_idx
  description: A length $p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.
returns: A $p \times p$ vector with $\mat{F}^{-1}_t$.
---


Update the precision of the forcast error.
Unlike `ssm_update_Finv`, this allows for missing values in `\vec{y}_{t}`.
If $y_{k,t}$ is missing, then $F^{-1}_{i,j,t} = 0$ for any $i = k$ or $j = k$,
otherwise it is the same as $\mat{F}^{-1}$ calculated after removing missing values.


This is the inverse of $\mat{F}_t$.

*/

matrix ssm_update_Finv_miss(matrix P, matrix Z, matrix H,
                                   int p_t, int[] y_idx) {
  matrix[rows(H), cols(H)] Finv;
  int p;
  int m;
  p = rows(H);
  m = cols(Z);
  if (p_t < p) {
    Finv = rep_matrix(0., p, p);
    if (p_t > 0) {
      matrix[p_t, m] Z_star;
      matrix[p_t, p_t] H_star;
      int idx[p_t];
      idx = y_idx[1:p_t];
      Z_star = Z[idx, :];
      H_star = H[idx, idx];
      Finv[idx, idx] = ssm_update_Finv(P, Z_star, H_star);
    }
  } else {
    Finv = ssm_update_Finv(P, Z, H);
  }
  return Finv;
}

/**
---
function: ssm_update_K
args:
- name: P
  description: An $m \times m$ vector with the variance of the predicted state, $P_t$.
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- name: T
  description: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- name: Finv
  description: A $p \times p$ matrix
returns: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
---

Update the Kalman gain, $\mat{K}_t$.


The Kalman gain is
$$
\mat{K}_t = \mat{T}_t \mat{P}_t \mat{Z}_t' \mat{F}^{-1}_t .
$$

*/

matrix ssm_update_K(matrix P, matrix Z, matrix T, matrix Finv) {
  matrix[cols(Z), rows(Z)] K;
  K = T * P * Z' * Finv;
  return K;
}


/**
---
function: ssm_update_L
args:
- name: Z
  description: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$
- name: T
  description: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- name: K
  description: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
returns: An $m \times m$ matrix, $\mat{L}_t$.
---

Update $L_t$


$$
\mat{L}_t = \mat{T}_t - \mat{K}_t \mat{Z}_t .
$$

*/

matrix ssm_update_L(matrix Z, matrix T, matrix K) {
  matrix[rows(T), cols(T)] L;
  L = T - K * Z;
  return L;
}


/**
---
function: ssm_update_loglik
args:
- name: v
  description: A $p \times 1$ matrix with the forecast error, $\vec{v}_t$.
- name: Finv
  description: A $p \times p$ matrix with variance of the forecast error, $\mat{F}^{-1}_t$.
returns: The log-likelihood
---


Calculate the log-likelihood of a single observation in a State-space model


The log-likehood of a single observation in a state-space model is
$$
\ell_t = - \frac{1}{2} p \log(2 \pi) - \frac{1}{2} \left(\log|\mat{F}_t| + \vec{v}_t' \mat{F}^{-1}_t \vec{v}_t  \right)
$$
*/

real ssm_update_loglik(vector v, matrix Finv) {
  real ll;
  int p;
  p = num_elements(v);
  // det(A^{-1}) = 1 / det(A) -> log det(A^{-1}) = - log det(A)
  ll = (- 0.5 *
        (p * log(2 * pi())
         - log_determinant(Finv)
         + quad_form_sym(Finv, v)
       ));
  return ll;
}

/**
---
function: ssm_update_loglik_miss
args:
- name: v
  description: A $p \times 1$ matrix with the forecast error, $\vec{v}_t$.
- name: Finv
  description: A $p \times p$ matrix with variance of the forecast error, $\mat{F}^{-1}_t$.
- name: p_t
  description: The number of non-missing elements in the observation vector, $\vec{y}_t$.
- name: y_idx
  description: A length $p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.
returns: A $p \times p$ vector with $\mat{F}^{-1}_t$.
---


Calculate the log-likelihood of a single observation in a State-space model

Unlike `ssm_update_loglik`, this allows for missing values.


*/

real ssm_update_loglik_miss(vector v, matrix Finv, int p_t, int[] y_idx) {
  real ll;
  int p;
  p = num_elements(v);
  if (p_t == 0) {
    ll = 0.;
  } else if (p_t == p) {
    ll = ssm_update_loglik(v, Finv);
  } else {
    int idx[p_t];
    matrix[p_t, p_t] Finv_star;
    vector[p_t] v_star;
    idx = y_idx[1:p_t];
    Finv_star = Finv[idx, idx];
    v_star = v[idx];
    ll = ssm_update_loglik(v_star, Finv_star);
  }
  return ll;
}


/**
`r render_section("Filtering")`

*/
/**
---
function: ssm_filter_idx
args:
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $6 \times 3$ integer array containing the indexes of the return values of the Kalman filter.
---

Indexes of the return values of the Kalman filter functions:
`ssm_filter`.


`ssm_filter_idx` returns a $6 \times 3$ integer array with the
(length, start index, stop index) of ($\ell_t$, $\vec{v}$, $\vec{F}^-1$, $\mat{K}$, $\vec{a}$, $\mat{P}$).

value            length                 start                               stop
---------------- ---------------------- ----------------------------------- ----------------------------------------------
$\ell_t$         $1$                    $1$                                 $1$
$\vec{v}$        $p$                    $2$                                 $1 + p$
$\mat{F}^{-1}$   $p (p + 1) / 2$        $2 + p$                             $1 + p + p (p + 1) / 2$
$\mat{K}$        $mp$                   $2 + p + p (p + 1) / 2$             $1 + p + p (p + 1) / 2 + mp$
$\vec{a}_t$      $m$                    $2 + p + p (p + 1) / 2 + mp$        $1 + p + p (p + 1) / 2 + mp + m$
$\mat{P}^t$      $m (m + 1) / 2$        $2 + p + p (p + 1) / 2 + mp + m$    $1 + p + p (p + 1) / 2 + mp + m (m + 1) / 2$


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

/**
---
function: ssm_filter_size
args:
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: The number of elements in the vector.
---

Number of elements in vector containing filter results


*/

int ssm_filter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  sz = idx[6, 3];
  return sz;
}

/**
---
function: ssm_filter_get_loglik
args:
- name: x
  description: A vector with results from `ssm_filter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: The log-likelihood $\ell_t$
---

Get the log-likehood from the results of `ssm_filter`.


*/

real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y = x[1];
  return y;
}

/**
---
function: ssm_filter_get_v
args:
- name: x
  description: vector with results from `ssm_filter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $p \times 1$ vector with the forecast error, $\vec{v}_t$.
---

Get the forecast error from the results of `ssm_filter`.


*/

vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[2, 2], idx[2, 1]);
  return y;
}

/**
---
function: ssm_filter_get_Finv
args:
- name: x
  description: vector with results from `ssm_filter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $p \times p$ matrix with the forecast precision, $\mat{F}^{-1}_t$.
---

Get the forecast precision from the results of `ssm_filter`.


*/

matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[3, 2], idx[3, 1]), p);
  return y;
}

/**
---
function: ssm_filter_get_K
args:
- name: x
  description: vector with results from `ssm_filter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: A $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.
---

Get the Kalman gain from the results of `ssm_filter`.


*/

matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[4, 2], idx[4, 1]), m, p);
  return y;
}

/**
---
function: ssm_filter_get_a
args:
- name: x
  description: vector with results from `ssm_filter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: An $m \times 1$ vector with the expected value of the predicted state, $\E(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \vec{a}_t$.
---

Get the expected value of the predicted state from the results of `ssm_filter`.


*/

vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[5, 2], idx[5, 1]);
  return y;
}

/**
---
function: ssm_filter_get_P
args:
- name: x
  description: vector with results from `ssm_filter`.
- name: m
  description: The number of states
- name: p
  description: The size of the observation vector $\vec{y}_t$.
returns: An $m \times m$ matrix with the variance of the predicted state, $\Var(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \mat{P}_t$.
---

Get the variance of the predicted state from the results of `ssm_filter`.


*/

matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[6, 2], idx[6, 1]), m);
  return y;
}

/**
---
function: ssm_filter
args:
- name: y
  description: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
returns: Array of size $n$ of $(1 + p + p (p + 1) / 2 + mp + m + m (m + 1) / 2) \times 1$ vectors in the format described in `ssm_filter_idx`.
---


Kalman filter


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

The results of Kalman filter for a given are returned as a $1 + p + p (p + 1) / 2 + m p + m (m + 1) / 2$ vector for each time period, where
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
  m = dims(Z)[3]; // number of states
  p = dims(Z)[2]; // obs size
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
    RQR = quad_form_sym(Q_t, R_t ');
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
        if (size(R) > 1 || size(Q) > 1) {
          RQR = quad_form_sym(Q_t, R_t ');
        }
      }
      // updating
      v = ssm_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_update_Finv(P, Z_t, H_t);
      K = ssm_update_K(P, Z_t, T_t, Finv);
      ll = ssm_update_loglik(v, Finv);
      // saving
      res[t, 1] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = symmat_to_vector(Finv);
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
  }
  return res;
}

/**
---
function: ssm_filter_miss
args:
- name: y
  description: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
- name: p_t
  description: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- name: y_idx
  description: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.
returns: A $p \times p$ vector with $\mat{F}^{-1}_t$.
---

The function `ssm_filter_miss` runs a Kalman filter like `ssm_filter`, except that it allows for missing values.

In the Kalman filter with missing values, the observation equation is,
$$
\begin{aligned}[t]
\vec{y}^*_{t} &= \vec{d}^*_{t} + \mat{Z}^*_{t} \vec{\alpha}_t + \vec{\varepsilon}^*_t \\
\vec{\varepsilon}^*_t &\sim N(\vec{0}, \mat{H}^*_t)
\end{aligned}
$$
where $\vec{y}^*_{t} = \mat{W}_t \vec{y}_t$, $\vec{d}^*_t = \mat{W}_t \vec{d}_t$,
$\mat{Z}^*_t = \mat{W}_t \mat{Z}_t$, $\mat{H}^*_t = \mat{W}_t \mat{H}_t \mat{W}_t\T$,
where $\mat{W}_t$ selects the non-missing rows of $\vec{y}_t$.

If all observations

If $y_{t,j}$ is missing, then

- $v_{t,j} = 0$
- $F^{-1}_{t,.,j} = \vec{0}_{p}$ and F^{-1}_{t,j,.} = \vec{0}_{p}$ if either $i
- $K_{t,.,j} = \vec{0}_{m}$

*/

vector[] ssm_filter_miss(vector[] y,
                          vector[] d, matrix[] Z, matrix[] H,
                          vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                          vector a1, matrix P1, int[] p_t, int[,] y_idx) {

  // returned data
  vector[ssm_filter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
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
    RQR = quad_form_sym(Q_t, R_t ');
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
        if (size(R) > 1 || size(Q) > 1) {
          RQR = quad_form_sym(Q_t, R_t ');
        }
      }
      // updating
      v = ssm_update_v_miss(y[t], a, d_t, Z_t, p_t[t], y_idx[t]);
      Finv = ssm_update_Finv_miss(P, Z_t, H_t, p_t[t], y_idx[t]);
      K = ssm_update_K(P, Z_t, T_t, Finv);
      ll = ssm_update_loglik_miss(v, Finv, p_t[t], y_idx[t]);
      // saving
      res[t, 1] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = symmat_to_vector(Finv);
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);
      // predict a_{t + 1}, P_{t + 1}
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
  }
  return res;
}

/**

`r render_section("Log-likelihood")`

*/
/**
---
function: ssm_lpdf
args:
- name: y
  description: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
returns: The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.
---

Log-likelihood of a Linear Gaussian State Space Model


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
  m = dims(Z)[3];
  p = dims(Z)[2];
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
    RQR = quad_form_sym(Q_t, R_t ');

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
          if (size(R) > 1 || size(Q) > 1) {
            RQR = quad_form_sym(Q_t, R_t ');
          }
        }
      }
      v = ssm_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_update_Finv(P, Z_t, H_t);
      K = ssm_update_K(P, T_t, Z_t, Finv);
      ll_obs[t] = ssm_update_loglik(v, Finv);
      // don't save a, P for last iteration
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}

/**
---
function: ssm_miss_lpdf
args:
- name: y
  description: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
- name: p_t
  description: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- name: y_idx
  description: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.
returns: The log-likelihood $p(\vec{y}_{1:n} | \vec{d}_{1:n}, \mat{Z}_{1:n}, \mat{H}_{1:n}, \vec{c}_{1:n}, \mat{T}_{1:n}, \mat{R}_{1:n}, \mat{Q}_{1:n})$.
---

*/
real ssm_miss_lpdf(vector[] y,
                   vector[] d, matrix[] Z, matrix[] H,
                   vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                   vector a1, matrix P1, int[] p_t, int[,] y_idx) {
  real ll;
  int n;
  int m;
  int p;
  int q;
  n = size(y); // number of obs
  m = dims(Z)[3];
  p = dims(Z)[2];
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
    RQR = quad_form_sym(Q_t, R_t ');

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
          if (size(R) > 1 || size(Q) > 1) {
            RQR = quad_form_sym(Q_t, R_t ');
          }
        }
      }
      v = ssm_update_v_miss(y[t], a, d_t, Z_t, p_t[t], y_idx[t]);
      Finv = ssm_update_Finv_miss(P, Z_t, H_t, p_t[t], y_idx[t]);
      K = ssm_update_K(P, Z_t, T_t, Finv);
      ll_obs[t] = ssm_update_loglik_miss(v, Finv, p_t[t], y_idx[t]);
      // don't save a, P for last iteration
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}


/**

`r render_section("Time-Invariant Kalman Filter")`

*/
/**
---
function: matrix_diff
args:
- name: A
  description: An $m \times n$ matrix.
- name: B
  description: An $m \times n$ matrix.
returns: If converged, then 1, else 0.
---

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

/**
---
function: ssm_constant_lpdf
args:
- name: y
  description: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
returns: The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.
---

Log-likelihood of a Time-Invariant Linear Gaussian State Space Model


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

    RQR = quad_form_sym(Q, R ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      v = ssm_update_v(y[t], a, d, Z);
      if (converged < 1) {
        Finv = ssm_update_Finv(P, Z, H);
        K = ssm_update_K(P, Z, T, Finv);
      }
      ll_obs[t] = ssm_update_loglik(v, Finv);
      // don't save a, P for last iteration
      if (t < n) {
        a = ssm_update_a(a, c, T, v, K);
        // check for convergence
        // should only check for convergence if there are no missing values
        if (converged < 1) {
          P_old = P;
          P = ssm_update_P(P, Z, T, RQR, K);
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

`r render_section("Common Smoother Functions")`

*/

/**
---
function: ssm_smooth_states_mean
args:
- name: filter
  description: The results of `ssm_filter`
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: c
  description: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- name: T
  description: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- name: R
  description: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- name: Q
  description: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
returns: An array of size $n$ of $m \times 1$ vectors containing $\hat{\vec{\alpha}}_t$.
---


The state smoother calculates $\hat{\vec{\alpha}}_t = \E(\vec{\alpha}_t | \vec{y}_{1:n})$.
$$
\hat{\vec{\alpha}}_{t + 1} = \mat{T}_t \hat{\vec{\alpha}}_{t} + \mat{R}_t \mat{Q}_t \mat{R}'_t \vec{r}_t ,
$$
where $r_t$ is calcualted from the state disturbance smoother.
The smoother is initialized at $t = 1$ with $\hat{\vec{\alpha}}_t = \vec{a}_1 + \mat{P}_1 \vec{r}_0$.

For  `Z`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]

*/

vector[] ssm_smooth_states_mean(vector[] filter,
                              matrix[] Z,
                              vector[] c, matrix[] T, matrix[] R, matrix[] Q) {
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
    vector[m] r[n + 1];
    vector[p] u;
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
    matrix[m, q] R_t;
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
      RQR = quad_form_sym(Q[1], R[1]');
    }
    // find smoothed state disturbances
    // Since I don't need to calculate the
    // variances of the smoothed disturbances,
    // I reimplement the state distrurbance smoother here
    // removing extraneous parts.
    // r goes from t = n, ..., 1, 0.
    // r_n
    r[n + 1] = rep_vector(0.0, m);
    for (s in 0:(n - 1)) {
      int t;
      // move backwards in time
      t = n - s;
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
      // u_{t - 1}
      u = Finv * v - K ' * r[t + 1];
      // r_{t - 1}
      r[t] = Z_t ' * u + T_t ' * r[t + 1];
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
        RQR = quad_form_sym(Q_t, R_t');
      }
      // `r[t + 1]` = $r_{t}$
      // alpha_{t + 1} = c_t + T_t * \alpha_t + R_t Q_t R'_t r_t
      alpha[t + 1] = c_t + T_t * alpha[t] + RQR * r[t + 1];
    }
  }
  return alpha;
}


/**
`r render_section("Simulators and Smoothing Simulators")`

*/
/**
---
function: ssm_sim_idx
args:
- name: m
  description: The number of states
- name: p
  description: The length of the observation vector
- name: q
  description: The number of state disturbances
returns: A 4 x 3 array of integers
---


Indexes of each component of `ssm_sim_rng` results.


The returned array has columns (length, start location, and end location)
for rows: $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\varepsilon}_t$, and $\vec{\eta}_t$ in the results of `ssm_sim_rng`.

*/

int[,] ssm_sim_idx(int m, int p, int q) {
  int sz[2, 3];
  // y
  sz[1, 1] = p;
  // a
  sz[2, 1] = m;
  // Fill in start and stop points
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  sz[2, 2] = sz[2 - 1, 3] + 1;
  sz[2, 3] = sz[2, 2] + sz[2, 1] - 1;
  return sz;
}

/**
---
function: ssm_sim_size
args:
- name: m
  description: The number of states
- name: p
  description: The length of the observation vector
- name: q
  description: The number of state disturbances
returns: The number of elements
---


The number of elements in vectors returned by `ssm_sim_rng` results.


*/

int ssm_sim_size(int m, int p, int q) {
  int sz;
  sz = ssm_sim_idx(m, p, q)[2, 3];
  return sz;
}

/**
---
function: ssm_sim_get_y
args:
- name: x
  description: vector of results from `ssm_sim_rng`
- name: m
  description: The number of states
- name: p
  description: The length of the observation vector
- name: q
  description: The number of state disturbances
returns: vector A $p \times 1$ vector with $\vec{y}_t$.
---

Extract $\vec{y}_t$ from vectors returned by `ssm_sim_rng`.

*/

vector ssm_sim_get_y(vector x, int m, int p, int q) {
  vector[p] y;
  int idx[2, 3];
  idx = ssm_sim_idx(m, p, q);
  y = x[idx[1, 2]:idx[1, 3]];
  return y;
}

/**
---
function: ssm_sim_get_a
args:
- name: x
  description: vector of results from `ssm_sim_rng`
- name: m
  description: The number of states
- name: p
  description: The length of the observation vector
- name: q
  description: The number of state disturbances
returns: A $m \times 1$ vector with $\vec{\alpha}_t$.
---

Extract $\vec{\alpha}_t$ from vectors returne by `ssm_sim_rng`.

*/

vector ssm_sim_get_a(vector x, int m, int p, int q) {
  vector[m] a;
  int idx[2, 3];
  idx = ssm_sim_idx(m, p, q);
  a = x[idx[2, 2]:idx[2, 3]];
  return a;
}


/**
---
function: ssm_sim_rng
args:
- name: n
  description: Number of time observations to draw, $t = 1, \dots, n$.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
returns: of size $n$ of vectors with Draw $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\eta}_t$ and $\vec{\varepsilon}_t$. See the description.
---

Simulate from a Linear Gaussian State Space model.

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

The returned vectors are of length $p + m$, in the format,
$$
(\vec{y}_t', \vec{\alpha}_t') .
$$
Use the functions `ssm_sim_get_y`, `ssm_sim_get_a` to extract values from the vector.

element         length         start         end
--------------- -------------- ------------- -----------
$y_t$           $p$            $1$           $p$
$\alpha$_t      $m$            $p + 1$       $p + m$

*/

vector[] ssm_sim_rng(int n,
                    vector[] d, matrix[] Z, matrix[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {
  vector[ssm_sim_size(dims(Z)[3], dims(Z)[2], dims(Q)[2])] ret[n];
  int p;
  int m;
  int q;
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    // system matrices for current iteration
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    matrix[p, p] HL;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[q, q] QL;
    // outputs
    vector[p] y;
    vector[p] eps;
    vector[m] a;
    vector[q] eta;
    // constants
    vector[p] zero_p;
    vector[q] zero_q;
    vector[m] zero_m;
    int idx[2, 3];

    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    HL = cholesky_decompose2(H_t);
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    QL = cholesky_decompose2(Q_t);

    idx = ssm_sim_idx(m, p, q);
    zero_p = rep_vector(0.0, p);
    zero_q = rep_vector(0.0, q);
    zero_m = rep_vector(0.0, m);
    a = multi_normal2_rng(a1, P1);
    for (t in 1:n) {
      // save alpha
      ret[t, idx[2, 2]:idx[2, 3]] = a;
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
          HL = cholesky_decompose2(H_t);
        }
      }
      // draw forecast error and observed value
      eps = multi_normal_cholesky2_rng(zero_p, HL);
      y = d_t + Z_t * a + eps;
      // save
      ret[t, idx[1, 2]:idx[1, 3]] = y;
      // calculate eta_{t} and a_{t + 1}
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
          QL = cholesky_decompose2(Q_t);
        }
        eta = multi_normal_cholesky2_rng(zero_q, QL);
        a = c_t + T_t * a + R_t * eta;
      }
    }
  }
  return ret;
}

/**

`r render_section("Simulation Smoothers")`

*/
/**
---
function: ssm_simsmo_states_rng
args:
- name: filter
  description: A length $n$ array with results from `ssm_filter`.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
returns: Array of size $n$ of $m \times 1$ vectors containing a single draw from $(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.
---

State simulation smoother

Draw samples from the posterior distribution of the states,
$\tilde{\vec{\alpha}}_{1:n} \sim p(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.


For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].

*/

vector[] ssm_simsmo_states_rng(vector[] filter,
                               vector[] d, matrix[] Z, matrix[] H,
                               vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                               vector a1, matrix P1) {
    vector[dims(Z)[3]] draws[size(filter)];
    int n;
    int p;
    int m;
    int q;
    n = size(filter);
    m = dims(Z)[3];
    p = dims(Z)[2];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter_plus[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[p] y[n];
      vector[m] alpha_hat_plus[n];
      vector[m] alpha_hat[n];
      // Smooth states
      alpha_hat = ssm_smooth_states_mean(filter, Z, c, T, R, Q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter with simulated y's
      filter_plus = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      alpha_hat_plus = ssm_smooth_states_mean(filter_plus, Z, c, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha_hat[i]);
      }
    }
    return draws;
}


/**
---
function: ssm_simsmo_states_miss_rng
args:
- name: filter
  description: A length $n$ array with results from `ssm_filter`.
- name: d
  description: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- name: Z
  description: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- name: H
  description: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
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
- name: p_t
  description: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- name: y_idx
  description: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.

returns: Array of size $n$ of $m \times 1$ vectors containing a single draw from $(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.
---

State simulation smoother, as in `ssm_simsmo_states_rng`, allowing for missing values.

*/

vector[] ssm_simsmo_states_miss_rng(vector[] filter,
                                    vector[] d, matrix[] Z, matrix[] H,
                                    vector[] c, matrix[] T,
                                    matrix[] R, matrix[] Q,
                                    vector a1, matrix P1,
                                    int[] p_t, int[,] y_idx) {
    vector[dims(Z)[3]] draws[size(filter)];
    int n;
    int p;
    int m;
    int q;
    n = size(filter);
    m = dims(Z)[3];
    p = dims(Z)[2];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter_plus[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[p] y[n];
      vector[m] alpha_hat_plus[n];
      vector[m] alpha_hat[n];
      // Smooth states
      alpha_hat = ssm_smooth_states_mean(filter, Z, c, T, R, Q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter with simulated y's
      filter_plus = ssm_filter_miss(y, d, Z, H, c, T, R, Q,
                                    a1, P1, p_t, y_idx);
      // mean correct epsilon samples
      alpha_hat_plus = ssm_smooth_states_mean(filter_plus, Z, c, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha_hat[i]);
      }
    }
    return draws;
}


/**

`r render_section("Stationary")`

*/
/**
---
function: pacf_to_acf
args:
- name: x
  description: A vector of coefficients of a partial autocorrelation function
returns: A vector of coefficients of an Autocorrelation function
---

Partial Autocorrelations to Autocorrelations

*/

// from R function partrans in arima.c
// https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/arima.c#L439
vector pacf_to_acf(vector x) {
  vector[num_elements(x)] x_new;
  vector[num_elements(x)] work;
  real a;
  int p;
  p = num_elements(x);
  work = x;
  x_new = x;
  if (p > 1) {
    for (j in 2:p) {
      a = x_new[j];
      for (k in 1:(j - 1)) {
        work[k] = work[k] - a * x_new[j - k];
      }
      for (k in 1:j) {
        x_new[k] = work[k];
      }
    }
  }
  return x_new;
}

/**
---
function: constrain_stationary
args:
- name: x
  description: An unconstrained vector in $(-\infty, \infty)$
returns: A vector of coefficients for a stationary AR or inverible MA process.
---

Constrain vector of coefficients to the stationary and intertible region for AR or MA functions.


See @Jones1980a, @Jones1987a, @Monahan1984a, @AnsleyKohn1986a, and the functions
`tools.constrain_stationary_univariate` and `tools.unconstraine_stationary_univariate` in
[statsmodels.tsa.statespace](http://www.statsmodels.org/dev/statespace.html#statespace-tools).


1. Each $\alpha_j \in (-\infty, \infty)$ is transformed to a partial correlation within $(-1, 1)$,
   with $\rho_j = \tanh(\alpha_j)$.
2. Then the partial correlations are converted to autocorrelation coefficients
   using the Durbin-Levinson recursions:
   $$
   $$

The transformation is reversed to take autocorrelation coefficients to
an unconstrained $R^p$ space.

1. Autocorrelation coefficients are transformed to partial autocorrelation coefficients,
    by running the Durbin-Levinson recursions in reverse.
2. Transform each partial autocorrelation to go from $(-1, 1) \to (-\infty, \infty)$ using,
    $\alpha_j = \atanh(\rho_j)$.



*/

vector constrain_stationary(vector x) {
  vector[num_elements(x)] r;
  int n;
  n = num_elements(x);
  // transform (-Inf, Inf) to (-1, 1)
  for (i in 1:n) {
    r[i] = tanh(x[i]);
  }
  // Transform PACF to ACF
  return pacf_to_acf(r);
}



/**
---
function: acf_to_pacf
args:
- name: x
  description: Coeffcients of an autocorrelation function.
returns: A vector of coefficients of the corresponding partial autocorrelation function.
---

Convert coefficients of an autocorrelation function to partial autocorrelations.

*/

// from R function invpartrans in arima.c
// https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/arima.c#L525
vector acf_to_pacf(vector x) {
  vector[num_elements(x)] x_new;
  vector[num_elements(x)] work;
  real a;
  int p;
  p = num_elements(x);
  work = x;
  x_new = x;
  // Run the Durbin-Levinson recursions backwards
  if (p > 1) {
    for(i in 0:(p - 2)) {
      int j;
      j = p - i;
      a = x_new[j];
      for(k in 1:(j - 1)) {
        work[k]  = (x_new[k] + a * x_new[j - k]) / (1 - pow(a, 2));
      }
      for (k in 1:j) {
        x_new[k] = work[k];
      }
    }
  }
  return x_new;
}

/**
---
function: unconstrain_stationary
args:
- name: x
  description: Coeffcients of an autocorrelation function.
returns: Coefficients of the corresponding partial autocorrelation function.
---

Transform from stationary and invertible space to $(-\infty, \infty)$.

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
    z[i] = atanh(r[i]);
  }
  return z;
}

/**
---
function: kronecker_prod
args:
- name: A
  description: An $m \times n$ matrix
- name: B
  description: A $p \times q$ matrix
returns: An $mp \times nq$ matrix.
---

The Kronecker product of a $\mat{A}$ and $\mat{B}$ is
$$
\mat{A} \otimes \mat{B} =
\begin{bmatrix}
a_{11} \mat{B} \cdots a_{1n} \mat{B} \\
\vdots & \ddots & vdots \\
a_{m1} \mat{B} & \cdots & a_{mn} \mat{B}
\end{bmatrix} .
$$

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

/**
---
function: stationary_cov
args:
- name: T
  description: The $m \times m$ transition matrix
- name: RQR
  description: The $m \times m$ system covarariance matrix, $\mat{R} \mat{Q} \mat{R}\T$.
returns: An $m \times m$ matrix with the stationary covariance matrix.
---


Find the covariance of the stationary distribution of an ARMA model


When a system is stationary, the initial covariance $P_1$ satistfies,
$$
\mat{P}_1 = \mat{T} \mat{P}_1 \mat{T}\T + \mat{R} \mat{Q} \mat{R}\T
$$
This can be solved as,
$$
(\mat{T} \otimes \mat{T}) \VEC(\mat{P}_1) = \VEC(\mat{R} \mat{Q} \mat{R}\T)
$$
where $\VEC(P_1)$ and $\VEC(R R')$ are the stacked columns of $\mat{P}_1$ and $\mat{R} \mat{Q} \mat{R}\T$

Note that in the special case of ARIMA models, $\mat{Q} = \sigma^2$, so
this simplifies to,
$$
(\mat{T} \otimes \mat{T}) \VEC(\mat{Q}_0) = \VEC(\mat{R} \mat{R}\T),
$$
where $P_1 = \sigma^2 \mat{Q}_0$.

In a stationary distribution, the initial mean $\vec{a}_1 = \vec{c}$.

See @DurbinKoopman2012, Sec 5.6.2.

*/

matrix stationary_cov(matrix T, matrix RQR) {
  matrix[rows(T), cols(T)] P;
  int m;
  m = rows(T);
  // m = 1 is an easy case, so treat it separately
  // since it doesn't require inverting a matrix
  if (m == 1) {
    P[1, 1] = RQR[1, 1] / (1.0 - pow(T[1, 1], 2));
  } else {
    matrix[rows(T) * rows(T), rows(T) * rows(T)] TT;
    vector[rows(T) * rows(T)] RQR_vec;
    int m2;
    m2 = m * m;
    RQR_vec = to_vector(RQR);
    # I_{m^2} - T \otimes T
    TT = - kronecker_prod(T, T);
    for (i in 1:m2) {
      TT[i, i] = 1.0 + TT[i, i];
    }
    P = to_matrix_colwise(inverse(TT) * RQR_vec, m, m);
  }
  return P;
}
