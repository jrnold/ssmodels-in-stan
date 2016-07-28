
#  Stan Functions

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

##  Utility Functions




### to_symmetric_matrix

**Arguments**

- `x`: An $n \times n$ matrix


**returns** An $n \times n$ symmetric matrix, $0.5 (x + x')$.

Ensure a matrix is symmetric


```stan
matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}


```



### to_matrix_colwise

**Arguments**

- `v`: An $n \times m$ vector.
- `m`: Number of rows in the vector
- `n`: Number of columns in the vector


**returns** A $m \times n$ matrix containting the elements from `v`

Convert vector to a matrix (column-major).


```stan
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


```



### matrix_pow

**Arguments**

- `A`: The matrix to take the power of
- `n`: The order of the power


**returns** None

Calculate the power of a matrix, $\mat{A}^n$.


```stan
matrix matrix_pow(matrix A, int n);

matrix matrix_pow(matrix A, int n) {
  if (n == 0) {
    return diag_matrix(rep_vector(1., rows(A)));
  } else if (n == 1) {
    return A;
  } else if (n > 1) {
    # recurively this is n log n.
    if (n % 2 == 0) {
      return matrix_pow(A, n / 2) * matrix_pow(A, n / 2);
    } else {
      return A * matrix_pow(A, n - 1);
    }
  } else {
    return A;
  }
}


```



### symmat_size

**Arguments**

- `n`: The number of rows and columns in the matrix.


**returns** The number of unique elements

Calculate the number of unique elements in a symmetric matrix

The number of unique elements in an $m \times m$ matrix is
$(m \times (m + 1)) / 2$.


```stan

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


```



### find_symmat_dim

**Arguments**

- `n`: The number of unique elements in a symmetric matrix.


**returns** The dimension of the associated symmetric matrix.

Given vector with $n$ elements containing the $m (m + 1) / 2$ elements of a symmetric matrix,
return $m$.


```stan
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


```



### vector_to_symmat

**Arguments**

- `x`: The vector with the unique elements
- `n`: The dimensions of the returned matrix, $n \times n$.


**returns** matrix An $n \times n$ symmetric matrix.

Convert a vector to a symmetric matrix


```stan
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


```



### symmat_to_vector

**Arguments**

- `x`: An $n \times n$ matrix.


**returns** A $n (n + 1) / 2$ vector with the unique elements in $x$.

Convert an $n \times n$ symmetric matrix to a length $n (n + 1) / 2$ vector
containing its unique elements.


```stan

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


```



### rep_lower_triangular_matrix

**Arguments**

- `x`: Value used for the non-zero elements of the matrix.
- `m`: number of rows
- `n`: number of columns
- `diag`: If true, then include 1's on the diagonal.


**returns** An $m \times n$ lower triangular matrix

Fill in an lower triangular matrix.


```stan

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


```



### rep_upper_triangular_matrix

**Arguments**

- `x`: Value used for the non-zero elements of the matrix.
- `m`: number of rows
- `n`: number of columns
- `diag`: If true, then include 1's on the diagonal.


**returns** An $m \times n$ upper triangular matrix

Fill in an upper triangular matrix


```stan

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


```



### rep_diagonal_matrix

**Arguments**

- `x`: Value used for the non-zero elements of the matrix.
- `m`: number of rows
- `n`: number of columns
- `k`: Index of the diagonal


**returns** An $m \times n$ upper triangular matrix

Create a diagonal $m \times n$ matrix with values $x$ on the $k$-th diagonal.


```stan

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


```



### fill_matrix

**Arguments**

- `x`: A $p \times q$, $p \leq m$, $\q \leq n$ matrix
- `m`: Number of rows in the returned matrix
- `n`: Number of columns in the returned matrix
- `i`: Indices mapping the rows of $A$ to the rows in the output matrix
- `j`: Indices mapping the columns of $A$ to the columns of the output matrix
- `a`: The default value in the returned matrix


**returns** A $m \times n$ matrix

Given a $p \times q$ matrix $\mat{X}$, default value $a$, and indexes $\vec{I} = i_1, ..., i_p$,
and $\vec{J} = j_1, ...j_q$, return a $m \times n$ matrix where $m \geq p$, $n \geq q$, where
$$
Y_{k, l} =
\begin{cases}
X_{i, j} & \text{if $k = i$, $l = j$, for some $i \in \vec{I}$, $j \in \vec{J}$,} \
a & \text{otherwise} .
\end{cases}
$$


```stan

matrix fill_matrix(matrix x, int m, int n, int[] i, int[] j, real a) {
  matrix[m, n] ret;
  ret = rep_matrix(a, m, n);
  ret[i, j] = x;
  return ret;
}


```



### fill_vector

**Arguments**

- `x`: A $p \times q$, $p \leq m$, $\q \leq n$ matrix
- `n`: Number of elements in the returned vector
- `i`: Indices mapping the rows of $A$ to the rows in the output matrix
- `a`: The default value in the returned vector


**returns** A $n \times 1$ matrix

Given an $m \times 1$ vector $\vec{x}$, an integer $n \geq m$, a default value $a$,
and indexes $\vec{I} = i_1, ..., i_m \in 1:n$, return a $n \times 1$ vector where
y_{j} =
\begin{cases}
x_{i} & \text{if $j = i$ for some $i \in \vec{I}$,} \
a & \text{otherwise}
\end{cases} .
$$


```stan

vector fill_vector(vector x, int n, int[] i, real a) {
  vector[n] ret;
  ret = rep_vector(a, n);
  ret[i] = x;
  return ret;
}


```



### int_sum_true

**Arguments**

- `x`: An array of length $n$ of integers


**returns** An integer between 0 and $n$.

For an array of integers, return the indexes where it is greater than zero.


```stan

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


```



### int_sum_false

**Arguments**

- `x`: An array of length $n$ of integers


**returns** An integer between 0 and $n$.

For an array of integers, return the indexes where it is less than or equal to zero.


```stan

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



```



### mask_indexes

**Arguments**

- `x`: An array of length $n$ of integers
- `n`: The number of false values in `x`.


**returns** An array of integers with elements having values between 1 and $m$.

For an array of integers, `x`, return the indexes where
mask is not true (`x[i] <= 0`).
The primary use of this function is where `x` represents
indicators for missing values,  and it is used to extract
the indexes of non-missing values.


```stan

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



```



### select_indexes

**Arguments**

- `x`: An array of length $m$ of integers
- `n`: The number of true values in `x`.


**returns** An array of integers with elements having values between 1 and $m$.

For an array of integers, `x`, return the indexes where
the elements are true (`x[i] > 0`).
The primary use of this function is where `x` represents
indicators for non-missing values, and it is used to extract
the indexes of non-missing values.


```stan

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


```



### normal2_rng

**Arguments**

- `mu`: mean
- `sigma`: variance


**returns** A value drawn from the specified normal distribution.

Draw samples from a normal distribution with mean `mu` and scale `sigma`.
Unlike the built-in `normal_rng()`, this allows for `sigma = 0`.


```stan

real normal2_rng(real mu, real sigma) {
  real y;
  if (sigma <= 0) {
    y = mu;
  } else {
    y = normal_rng(mu, sigma);
  }
  return y;
}


```



### cholesky_decompose2

**Arguments**

- `A`: An $n \times n$ matrix


**returns** An $n \times n$ lower-triangular matrix

Calculate the Cholesky decomposition of a matrix. Unlike the built-in
function, this handles cases in which the matrix has 0's on the diagonal.


```stan

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



```



### multi_normal2_rng

**Arguments**

- `mu`: An $n \times 1$ vector of the means
- `Sigma`: An $n \times n$ lower triangular matrix with covariance matrix.


**returns** An $n \times 1$ vector drawn from the specified multivariate normal distribution.

Sample from a multivariate normal distribution.
Unlike the built-in `multi_normal_rng`,
this function will still draw samples for deterministic elements in the vector.


```stan

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


```



### multi_normal_cholesky2_rng

**Arguments**

- `mu`: An $n \times 1$ vector of the means
- `L`: An $n \times n$ lower triangular matrix with the Cholesky decomposition of the covariance matrix.


**returns** An $n \times 1$ vector drawn from the specified multivariate normal distribution.

Sample from a multivariate normal distribution, parameterized with the Cholesky
decomposition of the covariance matrix. Unlike the built-in `multi_normal_cholesky_rng`,
this function will still draw samples for deterministic elements in the vector.


```stan

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



```


##  Filtering

Functions used in filtering and log-likelihood calculations.



### ssm_update_a

**Arguments**

- `a`: An $m \times 1$ vector with the predicted state, $\vec{a}_t$.
- `c`: An $m \times 1$ vector with the system intercept, $\vec{c}_t$
- `T`: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- `v`: A $p \times 1$ vector with the forecast error, $\vec{v}_t$.
- `K`: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.


**returns** A $m \times 1$ vector with the predicted state at $t + 1$, $\vec{a}_{t + 1}$.

Update the expected value of the predicted state, $\vec{a}_{t + 1} = \E(\vec{\alpha}_{t + 1} | \vec{y}_{1:t})$,


The predicted state $\vec{a}_{t + 1}$ is,
$$
\vec{a}_{t + 1} = \mat{T}_t \vec{a}_t + \mat{K}_t \vec{v}_t + \vec{c}_t .
$$


```stan

vector ssm_update_a(vector a, vector c, matrix T, vector v, matrix K) {
  vector[num_elements(a)] a_new;
  a_new = T * a + K * v + c;
  return a_new;
}


```



### ssm_update_a_u1

**Arguments**

- `a`: An $m \times 1$ vector with the predicted state, $\vec{a}_{t, i}$.
- `v`: The forecast error, $v_{t, i}$.
- `K`: A $m \times 1$ vector  with the Kalman gain, $\vec{K}_{t,i}$, for element $i$ of $\vec{y}_t$.


**returns** An $m \times 1$ vector with the predicted state $\vec{a}_{t, i + 1}$.

Update $\vec{a}_{t,i + 1}$ from $\vec{a}_{t, i}$ in the unvariate filter.


```stan
vector ssm_update_a_u1(vector a, real v, vector K) {
  vector[num_elements(a)] a_new;
  a_new = a + K * v;
  return a_new;
}


```



### ssm_update_a_u2

**Arguments**

- `a`: An $m \times 1$ vector with the predicted state, $\vec{a}_{t, p + 1}$.
- `c`: An $m \times 1$ vector with the system intercept, $\vec{c}_t$
- `T`: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.


**returns** An $m \times 1$ vector with the predicted state at $t + 1$, $\vec{a}_{t + 1,1}$.

Update the predicted state to the next time period, $\vec{a}_{t + 1,1}$ from $\vec{a}_{t, p + 1}$.


```stan
vector ssm_update_a_u2(vector a, vector c, matrix T) {
  vector[num_elements(a)] a_new;
  a_new = T * a + c;
  return a_new;
}




```



### ssm_update_P

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `T`: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- `RQR`: A $m \times m$ matrix with the system covariance matrix, $\mat{R}_t \mat{Q}_t \mat{R}_t'$.
- `K`: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.


**returns** An $m \times m$ matrix with the variance of the state, $\vec{P}_{t + 1}$.

Update the variance of the state in $t + 1$, $\mat{P}_{t + 1} = \Var(\alpha_{t + 1} | \vec{y}_{1:t})$,

The predicted state variance $\mat{P}_{t + 1}$ is,
$$
\mat{P}_{t + 1} = \mat{T}_t \mat{P}_t (\mat{T}_t - \mat{K}_t \mat{Z}_t)' + \mat{R}_t \mat{Q}_t \mat{R}_t' .
$$


```stan

matrix ssm_update_P(matrix P, matrix Z, matrix T,
                           matrix RQR, matrix K) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(T * P * (T - K * Z)' + RQR);
  return P_new;
}


```



### ssm_update_P_u1

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Finv`: The forecast precision, $F_{t,i}^{-1}$.
- `K`: An $m \times p$ matrix with the Kalman gain, $\mat{K}_{t,i}$.


**returns** An $m \times m$ matrix with the variance of the state, $\vec{P}_{t, i + 1}$.

Update the variance of the state in a univariate filter after observing $y_{t,i}$, $\mat{P}_{t, i + 1} = \Var(\alpha_{t} | \vec{y}_{1:t - 1}, y_{t,i}, \dots, y_{t,i})$.


```stan
matrix ssm_update_P_u1(matrix P, real Finv, vector K) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(P -  tcrossprod(to_matrix(K)) / Finv);
  return P_new;
}


```



### ssm_update_P_u2

**Arguments**

- `P`: An $m \times m$ vector with the variance of the state, $\mat{P}_{t, p + 1}$.
- `T`: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- `RQR`: A $m \times m$ matrix with the system covariance matrix, $\mat{R}_t \mat{Q}_t \mat{R}_t'$.


**returns** An $m \times m$ matrix with the variance of the state at $t + 1$ after observing $\vec{y_t}$, $\mat{P}_{t + 1,1}$.

Update the variance of the state in $t + 1$ in a univariate filter after observing $y_{t}$, $\mat{P}_{t, 1} = \Var(\alpha_{t + 1} | \vec{y}_{1:t})$.


```stan
matrix ssm_update_P_u2(matrix P, matrix T, matrix RQR) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(quad_form(P, T') + RQR);
  return P_new;
}


```



### ssm_update_v

**Arguments**

- `y`: A $p \times 1$ vector of the observations, $\vec{y}_t$.
- `a`: A $m \times 1$ vector of the states, $\vec{a}_t$.
- `d`: A $p \times 1$ vector with the observation intercept, $\vec{d}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.


**returns** A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.

Update the forcast error, $\vec{v}_t = \vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}})$


The forecast error $\vec{v}_t$ is
$$
\vec{v}_t =\vec{y}_t - \mat{Z}_t \vec{a}_t - \vec{d}_t .
$$


```stan

vector ssm_update_v(vector y, vector a, vector d, matrix Z) {
  vector[num_elements(y)] v;
  v = y - Z * a - d;
  return v;
}



```



### ssm_update_v_u

**Arguments**

- `y`: An observation, $y_{t, j}$.
- `a`: A $m \times 1$ vector of the states, $\vec{a}_t$.
- `d`: The observation intercept, $d_{t,j}$.
- `Z`: A $1 \times m$ vector from the design matrix, $\vec{Z}_{t, j}$.


**returns** The forecast errors, $v_{t,j}$.

Update the forcast error in univariate filtering.


```stan
real ssm_update_v_u(real y, vector a, real d, row_vector Z) {
  real v;
  v = y - dot_product(Z, a) - d;
  return v;
}


```



### ssm_update_v_miss

**Arguments**

- `y`: A $p \times 1$ vector of the observations, $\vec{y}_t$.
- `a`: A $m \times 1$ vector of the states, $\vec{a}_t$.
- `d`: A $p \times 1$ vector with the observation intercept, $\vec{d}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `p_t`: The number of non-missing elements of the observation vector $\vec{y}_t$.
- `y_idx`: A length $p$ array of integers indexes with the indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.

Update the forcast error, but unlike `ssm_update_v`, allow for missing
values in $\vec{y}_t$.


The elements of the forecast error $\vec{v}_t$ is
$$
\vec{v}_t =
\begin{cases}
  y_{j,t} - \vec{Z}_{j,.,t} \vec{a}_t - d_{j,t} & \text{if $y_{j,t} not missing.} \
  0 & \text{if $y_{j,t}$ is missing.}
\end{cases}
$$


```stan

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


```



### ssm_update_F

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `H`: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.


**returns** A $p \times p$ vector with $\mat{F}_t$.

Update the variance of the forcast error, $\mat{F}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))$


The variance of the forecast error $\mat{F}_t$ is
$$
\mat{F}_t = \mat{Z}_t \mat{P}_t \mat{Z}_t + \mat{H}_t .
$$


```stan

matrix ssm_update_F(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] F;
  F = to_symmetric_matrix(quad_form(P, Z') + H);
  return F;
}


```



### ssm_update_Finv

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `H`: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.


**returns** A $p \times p$ vector with $\mat{F}^{-1}_t$.

Update the precision of the forcast error, $\mat{F}^{-1}_t = \Var(\vec{y}_t - \E(\vec{y}_t | \vec{y_{1:(t - 1)}}))^{-1}$


This is the inverse of $\mat{F}_t$.


```stan

matrix ssm_update_Finv(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] Finv;
  // if can guarantee that F is spd, then take spd inverse.
  Finv = inverse_spd(to_symmetric_matrix(quad_form(P, Z') + H));
  // Finv = inverse(quad_form(P, Z') + H);
  return Finv;
}


```



### ssm_update_F_u

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times 1$ row vector from the design matrix, $\mat{Z}_{t,i}$.
- `H`: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.


**returns** The univariate forecast error, $f_{t,j}$.

Update the variance of the univariate forcast error, $f^{-1}_{t,i}$.


```stan
real ssm_update_F_u(matrix P, row_vector Z, real H) {
  return quad_form(P, Z') + H;
}



```



### ssm_update_Finv_u

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times 1$ row vector from the design matrix, $\mat{Z}_{t,j}$.
- `H`: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.


**returns** The univariate forecast precision, $f^{-1}_{t,j}$.

Update the precision of the univariate forcast error, $f^{-1}_{t,j}$.


```stan
real ssm_update_Finv_u(matrix P, row_vector Z, real H) {
  return 1. / (quad_form(P, Z') + H);
}


```



### ssm_update_Finv_miss

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `H`: A $p \times p$ matrix with the observation covariance matrix, $\mat{H}_t$.
- `p_t`: The number of non-missing elements in the observation vector, $\vec{y}_t$.
- `y_idx`: A length $p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** A $p \times p$ vector with $\mat{F}^{-1}_t$.

Update the precision of the forcast error.
Unlike `ssm_update_Finv`, this allows for missing values in `\vec{y}_{t}`.
If $y_{k,t}$ is missing, then $F^{-1}_{i,j,t} = 0$ for any $i = k$ or $j = k$,
otherwise it is the same as $\mat{F}^{-1}$ calculated after removing missing values.


This is the inverse of $\mat{F}_t$.


```stan

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


```



### ssm_update_K

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $P_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `T`: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- `Finv`: A $p \times p$ matrix


**returns** An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.

Update the Kalman gain, $\mat{K}_t$.


The Kalman gain is
$$
\mat{K}_t = \mat{T}_t \mat{P}_t \mat{Z}_t' \mat{F}^{-1}_t .
$$


```stan

matrix ssm_update_K(matrix P, matrix Z, matrix T, matrix Finv) {
  matrix[cols(Z), rows(Z)] K;
  K = T * P * Z' * Finv;
  return K;
}


```



### ssm_update_K_u

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $P_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `Finv`: $ matrix


**returns** An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.

Update the Kalman gain, $\mat{K}_t$, in univiariate filtering.


```stan

vector ssm_update_K_u(matrix P, row_vector Z, real Finv) {
  vector[num_elements(Z)] K;
  K = P * Z' * Finv;
  return K;
}



```



### ssm_update_L

**Arguments**

- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$
- `T`: An $m \times m$ matrix with the transition matrix, $\mat{T}_t$.
- `K`: An $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.


**returns** An $m \times m$ matrix, $\mat{L}_t$.

Update $L_t$


$$
\mat{L}_t = \mat{T}_t - \mat{K}_t \mat{Z}_t .
$$


```stan

matrix ssm_update_L(matrix Z, matrix T, matrix K) {
  matrix[rows(T), cols(T)] L;
  L = T - K * Z;
  return L;
}


```



### ssm_update_L_u

**Arguments**

- `Z`: A $1 \times m$ row vector from the design matrix, $\vec{Z}_{t,i}$
- `K`: An $m \times 1$ matrix with the Kalman gain, $\vec{K}_{t,i}$.


**returns** An $m \times m$ matrix, $\mat{L}_{t,i}$.

Update $L_t$ for univariate filtering,
$$
\mat{L}_t = \mat{I}_m - \vec{K}_{t,i}\vec{Z}_{t,i}
$$
See []@DurbinKoopman2012, p. 157]


```stan

matrix ssm_update_L_u(row_vector Z, vector K) {
  matrix[num_elements(Z), num_elements(Z)] L;
  int m;
  m = num_elements(Z);
  L = diag_matrix(rep_vector(1., m)) - K * Z;
  return L;
}


```



### ssm_update_loglik

**Arguments**

- `v`: A $p \times 1$ matrix with the forecast error, $\vec{v}_t$.
- `Finv`: A $p \times p$ matrix with variance of the forecast error, $\mat{F}^{-1}_t$.


**returns** The log-likelihood

Calculate the log-likelihood of a single observation in a State-space model


The log-likehood of a single observation in a state-space model is
$$
\ell_t = - \frac{1}{2} p \log(2 \pi) - \frac{1}{2} \left(\log|\mat{F}_t| + \vec{v}_t' \mat{F}^{-1}_t \vec{v}_t  \right)
$$


```stan

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


```



### ssm_update_loglik_miss

**Arguments**

- `v`: A $p \times 1$ matrix with the forecast error, $\vec{v}_t$.
- `Finv`: A $p \times p$ matrix with variance of the forecast error, $\mat{F}^{-1}_t$.
- `p_t`: The number of non-missing elements in the observation vector, $\vec{y}_t$.
- `y_idx`: A length $p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** A $p \times p$ vector with $\mat{F}^{-1}_t$.

Calculate the log-likelihood of a single observation in a State-space model

Unlike `ssm_update_loglik`, this allows for missing values.


```stan

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


```



### ssm_update_loglik_u

**Arguments**

- `v`: The forecast error, $v_{t,i}$.
- `Finv`: The forecast error, $f^{-1}_{t,i}$.


**returns** The log-likelihood

Calculate the log-likelihood of a single observation using univariate filtering.


```stan

real ssm_update_loglik_u(real v, real Finv) {
  real ll;
  // det(A^{-1}) = 1 / det(A) -> log det(A^{-1}) = - log det(A)
  ll = (- 0.5 * (
         log(2 * pi())
         - log(Finv)
         + Finv * pow(v, 2)
       ));
  return ll;
}



```

##  Filtering




### ssm_filter_idx

**Arguments**

- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** A $6 \times 3$ integer array containing the indexes of the return values of the Kalman filter.

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


```stan

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


```



### ssm_filter_size

**Arguments**

- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** The number of elements in the vector.

Number of elements in vector containing filter results


```stan

int ssm_filter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  sz = idx[6, 3];
  return sz;
}


```



### ssm_filter_get_loglik

**Arguments**

- `x`: A vector with results from `ssm_filter`.
- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** The log-likelihood $\ell_t$

Get the log-likehood from the results of `ssm_filter`.


```stan

real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y = x[1];
  return y;
}


```



### ssm_filter_get_v

**Arguments**

- `x`: vector with results from `ssm_filter`.
- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** A $p \times 1$ vector with the forecast error, $\vec{v}_t$.

Get the forecast error from the results of `ssm_filter`.


```stan

vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[2, 2], idx[2, 1]);
  return y;
}


```



### ssm_filter_get_Finv

**Arguments**

- `x`: vector with results from `ssm_filter`.
- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** A $p \times p$ matrix with the forecast precision, $\mat{F}^{-1}_t$.

Get the forecast precision from the results of `ssm_filter`.


```stan

matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[3, 2], idx[3, 1]), p);
  return y;
}


```



### ssm_filter_get_K

**Arguments**

- `x`: vector with results from `ssm_filter`.
- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** A $m \times p$ matrix with the Kalman gain, $\mat{K}_t$.

Get the Kalman gain from the results of `ssm_filter`.


```stan

matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[4, 2], idx[4, 1]), m, p);
  return y;
}


```



### ssm_filter_get_a

**Arguments**

- `x`: vector with results from `ssm_filter`.
- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** An $m \times 1$ vector with the expected value of the predicted state, $\E(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \vec{a}_t$.

Get the expected value of the predicted state from the results of `ssm_filter`.


```stan

vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[5, 2], idx[5, 1]);
  return y;
}


```



### ssm_filter_get_P

**Arguments**

- `x`: vector with results from `ssm_filter`.
- `m`: The number of states
- `p`: The size of the observation vector $\vec{y}_t$.


**returns** An $m \times m$ matrix with the variance of the predicted state, $\Var(\vec{\alpha}_t | \vec{y}_{1:(t-1)}) = \mat{P}_t$.

Get the variance of the predicted state from the results of `ssm_filter`.


```stan

matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[6, 2], idx[6, 1]), m);
  return y;
}


```



### ssm_filter

**Arguments**

- `y`: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** Array of size $n$ of $(1 + p + p (p + 1) / 2 + mp + m + m (m + 1) / 2) \times 1$ vectors in the format described in `ssm_filter_idx`.

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


```stan

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


```



### ssm_filter_miss

**Arguments**

- `y`: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
- `p_t`: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- `y_idx`: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** A $p \times p$ vector with $\mat{F}^{-1}_t$.

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


```stan

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


```



### ssm_filter_states_size

**Arguments**

- `m`: Number of states


**returns** The size of the vector

Length of the vectors returned by `ssm_filter_states`.


```stan

int ssm_filter_states_size(int m) {
  int sz;
  sz = m + symmat_size(m);
  return sz;
}


```



### ssm_filter_states_get_a

**Arguments**

- `x`: A vector returned by `ssm_filter_states`
- `m`: Number of states


**returns** An $m \times 1$ vector with the filtered expected value of the state, $\vec{a}_{t|t} = \E(\vec{\alpha}_t | \vec{y}_{1:t})$.

Extract $a_{t|t}$ from the results of `ssm_filter_states`


```stan

vector ssm_filter_states_get_a(vector x, int m) {
  vector[m] a;
  a = x[ :m];
  return a;
}


```



### ssm_filter_states_get_P

**Arguments**

- `x`: A vector returned by `ssm_filter_states`
- `m`: Number of states


**returns** An $m \times m$ matrix with the filtered variance of the state, $\mat{P}_{t|t} = \Var(\vec{\alpha}_t | \vec{y}_{1:t})$.

Extract $P_{t|t}$ from the results of `ssm_filter_states`


```stan

matrix ssm_filter_states_get_P(vector x, int m) {
  matrix[m, m] P;
  P = vector_to_symmat(x[(m + 1): ], m);
  return P;
}


```



### ssm_filter_states_update_a

**Arguments**

- `a`: An $m \times 1$ vector with the expected value of the predicted state, $\vec{a}_t$.
- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `v`: A $p \times 1$ vector with the forecast errors, $\vec{v}_t$
- `Finv`: A $p \times p$ matrix with the forecast prediction, $\mat{F}_{t}^{-1}$.


**returns** An $m \times 1$ matrix the expected value of the fitered states, $\E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t}$.

Calculate filtered state values [@DurbinKoopman2012, Sec 4.3.2],
$$
\E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t} = \mat{T}_t \vec{a}_t + \mat{K}_t \vec{v}_t .
$$


```stan

vector ssm_filter_states_update_a(vector a, matrix P, matrix Z,
                                  vector v, matrix Finv) {
  vector[num_elements(a)] aa;
  aa = a + P * Z ' * Finv * v;
  return aa;
}


```



### ssm_filter_states_update_P

**Arguments**

- `P`: An $m \times m$ vector with the variance of the predicted state, $\mat{P}_t$.
- `Z`: A $p \times m$ matrix with the design matrix, $\mat{Z}_t$.
- `Finv`: A $p \times p$ matrix with the forecast prediction, $\mat{F}_{t}^{-1}$.


**returns** An $m \times m$ matrix with variance of the filtered states, $\Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t}$.

Calculate filtered state variance values [@DurbinKoopman2012, Sec 4.3.2],
$$
\Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t} = \mat{P}_t - \mat{P}_t \mat{Z}_t' \mat{F}_t^{-1} \mat{Z}_t \mat{P}_t .
$$


```stan

matrix ssm_filter_states_update_P(matrix P, matrix Z, matrix Finv) {
  matrix[rows(P), cols(P)] PP;
  PP = to_symmetric_matrix(P - P * quad_form_sym(Finv, Z) * P);
  return PP;
}



```



### ssm_filter_states

**Arguments**

- `filter`: Results from `ssm_filter`
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.


**returns** of size $n$ of vectors.

Calculate filtered expected values and variances of the states

The filtering function `ssm_filter` returns the mean and variance of the predicted states,
$\vec{a}_t = \E(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$ and $\mat{P}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:(t -1)})$.


The vectors returned by `ssm_filter_states` are of length $m + m ^ 2$, with
$$
\vec{v}_t = (\vec{a}_{t|t}', \VEC(\vec{P}_{t|t})' )'
$$
Use the functions `ssm_filter_states_get_a` and `ssm_filter_states_get_P` to extract
elements from the results.

For `Z` the array can have a size of 1, if it is not time-varying, or a size of $n - 1$ if it is time varying.


```stan

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
      a = ssm_filter_states_update_a(a, P, Z_t, v, Finv);
      P = ssm_filter_states_update_P(P, Z_t, Finv);
      // saving
      res[t, :m] = a;
      res[t, (m + 1): ] = symmat_to_vector(P);
    }
  }
  return res;
}


```


##  Log-likelihood




### ssm_lpdf

**Arguments**

- `y`: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.

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


```stan

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


```



### ssm_miss_lpdf

**Arguments**

- `y`: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
- `p_t`: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- `y_idx`: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** The log-likelihood $p(\vec{y}_{1:n} | \vec{d}_{1:n}, \mat{Z}_{1:n}, \mat{H}_{1:n}, \vec{c}_{1:n}, \mat{T}_{1:n}, \mat{R}_{1:n}, \mat{Q}_{1:n})$.




```stan
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


```



### ssm_ufilter_lpdf

**Arguments**

- `y`: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.

Log-likelihood of a Linear Gaussian State Space Model calculated using
univariate filtering.


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


```stan

real ssm_ufilter_lpdf(vector[] y,
                        vector[] d, matrix[] Z, vector[] H,
                        vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                        vector a1, matrix P1) {
  // returned data
  real ll;
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
    real v;
    real Finv;
    vector[m] K;
    matrix[n, p] ll_obs;

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
      // update step
      for (i in 1:p) {
        Z_ti = row(Z_t, i);
        v = ssm_update_v_u(y[t, i], a, d_t[i], Z_ti);
        Finv = ssm_update_Finv_u(P, Z_ti, H_t[i]);
        K = ssm_update_K_u(P, Z_ti, Finv);
        ll[t, i] = ssm_update_loglik_u(v, Finv);
        a = ssm_update_a_u1(a, v, K);
        P = ssm_update_P_u1(P, Finv, K);
      }
      // predict step
      if (t < n) {
        a = ssm_update_a_u2(a, c_t, T_t);
        P = ssm_update_P_u2(P, T_t, RQR);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}

real ssm_ufilter_miss_lpdf(vector[] y,
                          vector[] d, matrix[] Z, vector[] H,
                          vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                          vector a1, matrix P1, int[,] miss) {
  // returned data
  real ll;
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
    real v;
    real Finv;
    vector[m] K;
    matrix[n, p] ll_obs;

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
      // update step
      for (i in 1:p) {
        if (miss[i, p]) {
          ll_obs[t, i] = 0.;
        } else {
          Z_ti = row(Z_t, i);
          v = ssm_update_v_u(y[t, i], a, d_t[i], Z_ti);
          Finv = ssm_update_Finv_u(P, Z_ti, H_t[i]);
          K = ssm_update_K_u(P, Z_ti, Finv);
          ll_obs[t, i] = ssm_update_loglik_u(v, Finv);
          a = ssm_update_a_u1(a, v, K);
          P = ssm_update_P_u1(P, Finv, K);
        }
      }
      // predict step
      if (t < n) {
        a = ssm_update_a_u2(a, c_t, T_t);
        P = ssm_update_P_u2(P, T_t, RQR);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}



```


##  Time-Invariant Kalman Filter




### matrix_diff

**Arguments**

- `A`: An $m \times n$ matrix.
- `B`: An $m \times n$ matrix.


**returns** If converged, then 1, else 0.

The difference between $A$ and $B$ is calculated as,
$$
d(A, B) = \max(A - B) / \max(A)
$$


```stan

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


```



### ssm_constant_lpdf

**Arguments**

- `y`: Observations, $\vec{y}_t$. An array of size $n$ of $p \times 1$ vectors.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** The log-likelihood, $p(\vec{y}_{1:n} | \vec{d}, \mat{Z}, \mat{H}, \vec{c}, \mat{T}, \mat{R}, \mat{Q})$, marginalized over the latent states.

Log-likelihood of a Time-Invariant Linear Gaussian State Space Model


Unlike `ssm_filter`, this function requires the system matrices (`d`, `Z`, `H`, `c`, `T`, `R`, `Q`)
to all be time invariant (constant).
When the state space model is time-invariant, then the Kalman recursion for $\mat{P}_t$ converges.
This function takes advantage of this feature and stops updating $\mat{P}_t$ after it converges
to a steady state.


```stan

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




```


##  Common Smoother Functions




### ssm_update_r

**Arguments**

- `r`: An $m \times 1$ vector with $\vec{r}_{t-1}$
- `Z`: A $p \times m$ vector with $\mat{Z}_{t}$
- `v`: A $p \times 1$ vector of the forecast errors, $\vec{v}_t$.
- `Finv`: A $p \times p$ matrix of the forecast precision, $\mat{F}^{-1}_t$.
- `L`: An $m \times m$ matrix with $\mat{L}_t$.


**returns** An $m \times 1$ vector with $\vec{r}_t$.

Update $\vec{r}_t$ in smoothing recursions

In smoothing recursions, the vector $\vec{r}_t$ is updated with,
$$
\vec{r}_{t - 1} = \mat{Z}' \mat{F}^{-1}_t \vec{v}_t + \mat{L}' \vec{r}_{t} .
$$

See [@DurbinKoopman2012, Sec 4.4.4, p. 91]


```stan

vector ssm_update_r(vector r, matrix Z, vector v, matrix Finv,
                           matrix L) {
  vector[num_elements(r)] r_new;
  r_new = Z ' * Finv * v + L ' * r;
  return r_new;
}


```


##  Common Smoother Functions




### ssm_update_r_u1

**Arguments**

- `r`: An $m \times 1$ vector with $\vec{r}_{t,i}$
- `Z`: A $1 \times m$ row vector of the design matix $\mat{Z}_{t,i}$
- `v`: The forecast error, $\vec{v}_{t,i}$.
- `Finv`: The $forecast precision, $\mat{F}^{-1}_{t,i}$.
- `L`: An $m \times m$ matrix with $\mat{L}_{t,i}$.


**returns** An $m \times 1$ vector with $\vec{r}_{t,i-1}$.

Update $\vec{r}_{t,i-1}$ from $\vec{r}_{t,i}$ in univariate smoothing recursions.

See [@KoopmanDurbin2012, p. 157]


```stan

vector ssm_update_r_u1(vector r, row_vector Z, real v, real Finv, matrix L) {
  vector[num_elements(r)] r_new;
  r_new = Z ' * Finv * v + L ' * r;
  return r_new;
}


```



### ssm_update_r_u2

**Arguments**

- `r`: An $m \times 1$ vector with $\vec{r}_{t,i}$
- `T`: A $m \times m$ row vector with the transition matrix $\mat{T}_{t-1}$


**returns** An $m \times 1$ vector with $\vec{r}_{t - 1,p}$.

Update $\vec{r}_{t-1,p}$ from $\vec{r}_{t,0}$ in univariate smoothing recursions.

See [@KoopmanDurbin2012, p. 157]


```stan

vector ssm_update_r_u2(vector r, matrix T) {
  vector[num_elements(r)] r_new;
  r_new = T ' * r;
  return r_new;
}


```



### ssm_update_N

**Arguments**

- `N`: An $m \times 1$ vector with $\vec{N}_{t-1}$
- `Z`: A $p \times m$ vector with $\mat{Z}_{t}$
- `Finv`: A $p \times p$ matrix of the forecast precision, $\mat{F}^{-1}_t$.
- `L`: An $m \times m$ matrix with $\mat{L}_t$.


**returns** An $m \times m$ matrix with $\vec{N}_t$.

Update $\mat{N}_t$ in smoothing recursions

In smoothing recursions, the matrix $\vec{N}_t$ is updated with,
$$
\mat{N}_{t - 1} = \mat{Z}_t' \mat{F}^{-1}_t \mat{Z}_t + \mat{L}_t' \mat{N}_t \mat{L}_t .
$$

See [@DurbinKoopman2012, Sec 4.4.4, p. 91]


```stan

matrix ssm_update_N(matrix N, matrix Z, matrix Finv, matrix L) {
  matrix[rows(N), cols(N)] N_new;
  # may not need this to_symmetric_matrix
  N_new = quad_form_sym(Finv, Z) + quad_form_sym(N, L);
  return N_new;
}


```



### ssm_update_N_u1

**Arguments**

- `N`: An $m \times 1$ vector with $\vec{N}_{t,i}$
- `Z`: A $1 \times m$ vector with $\mat{Z}_{t,i}$
- `Finv`: The forecast precision, $\mat{F}^{-1}_{t,i}$.
- `L`: An $m \times m$ matrix with $\mat{L}_t$.


**returns** An $m \times m$ matrix with $\vec{N}_t$.

Filter $\mat{N}_{t,i}$ after observing $y_{t,i}$ in univariate smoothing recursions,
$$
\mat{N}_{t,i-1} = \mat{Z}_{t,i}' \mat{F}^{-1}_{t,i} \mat{Z}_{t,i} + \mat{L}_{t,i}' \mat{N}_{t,i} \mat{L}_{t,i} .
$$

See [@DurbinKoopman2012, Eq. 6.15, p. 157]


```stan

matrix ssm_update_N_u1(matrix N, row_vector Z, real Finv, matrix L) {
  matrix[rows(N), cols(N)] N_new;
  # may not need this to_symmetric_matrix
  N_new = crossprod(to_matrix(Z)) * Finv + quad_form_sym(N, L);
  return N_new;
}


```



### ssm_update_N_u2

**Arguments**

- `N`: An $m \times 1$ vector with $\vec{N}_{t,0}$
- `T`: The $m \times m$ transition matrix $\mat{T}_{t - 1}$


**returns** An $m \times m$ matrix with $\vec{N}_{t-1,p}

Update smoothing variance from $t$ to $t - 1$, $\mat{N}_{t,0}$ to $\mat{N}_{t - 1, p}
$$
\mat{N}_{t,i-1} =  \mat{T}_{t-1}' \mat{N}_{t,0} \mat{T}_{t-1} .
$$

See [@DurbinKoopman2012, Eq. 6.15, p. 157]


```stan

matrix ssm_update_N_u2(matrix N, matrix T) {
  matrix[rows(N), cols(N)] N_new;
  # may not need this to_symmetric_matrix
  N_new = quad_form_sym(N, T);
  return N_new;
}



```



### ssm_smooth_state_size

**Arguments**

- `m`: The number of states.


**returns** The size of the vectors is $m + m (m + 1) / 2$.

The number of elements in vectors returned by `ssm_smooth_state`


```stan

int ssm_smooth_state_size(int m) {
  int sz;
  sz = m + symmat_size(m);
  return sz;
}


```



### ssm_smooth_state_get_mean

**Arguments**

- `x`: A vector returned by `ssm_smooth_state`
- `m`: The size of the state vector, $\vec{\alpha}_t$.


**returns** An $m \times 1$ vector with $\hat{\vec{\alpha}}_t$.

Extract $\hat{\vec{\alpha}}_t$ from vectors returned by `ssm_smooth_state`


```stan

vector ssm_smooth_state_get_mean(vector x, int m) {
  vector[m] alpha;
  alpha = x[ :m];
  return alpha;
}


```



### ssm_smooth_state_get_var

**Arguments**

- `x`: A vector returned by `ssm_smooth_state`
- `m`: The number of states


**returns** An $m \times m$ matrix with $\mat{V}_t$.

Extract $mat{V}_t$ from vectors returned by `ssm_smooth_state`


```stan

matrix ssm_smooth_state_get_var(vector x, int m) {
  matrix[m, m] V;
  V = vector_to_symmat(x[(m + 1): ], m);
  return V;
}



```



### ssm_smooth_state

**Arguments**

- `filter`: Results of `ssm_filter`
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.


**returns** An array of vectors constaining $\hat{\vec{\alpha}}_t$ and $\mat{V}_t = \Var(\vec{\alpha}_t | \vec{y}_{1:n})$.

The state smoother

This calculates the mean and variance of the states, $\vec{\alpha}_t$, given the entire sequence, $\vec{y}_{1:n}$.

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
$\hat{\vec{\alpha}}_t$                    $m$             $1$                   $m$
$\mat{V}_t$                                 $m (m + 1) / 2$ $m + 1$               $m + m (m + 1) / 2$


See @DurbinKoopman2012, Eq 4.44 and eq 4.69.


```stan

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
    // results
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
      L = ssm_update_L(Z_t, T_t, K);
      // r_{t - 1} and N_{t - 1}
      r = ssm_update_r(r, Z_t, v, Finv, L);
      N = ssm_update_N(N, Z_t, Finv, L);
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



```



### ssm_smooth_eps_size

**Arguments**

- `p`: The length of the observation vectors, $\vec{y}_t$.


**returns** The size of the vectors is $p + p (p + 1) / 2$.

The size of the vectors returned by `ssm_smooth_eps`


```stan

int ssm_smooth_eps_size(int p) {
  int sz;
  sz = p + symmat_size(p);
  return sz;
}


```



### ssm_smooth_eps_get_mean

**Arguments**

- `x`: vector from the results of `ssm_smooth_eps`.
- `p`: The length of the observation vectors, $\vec{y}_t$.


**returns** A $p \times 1$ vector with $\hat{\vec{\varepsilon}}_t$.

Extract $\hat{\vec{\varepsilon}}_t$ from vectors returned by `ssm_smooth_eps`


```stan

vector ssm_smooth_eps_get_mean(vector x, int p) {
  vector[p] eps;
  eps = x[ :p];
  return eps;
}


```



### ssm_smooth_eps_get_var

**Arguments**

- `x`: A vector returned by `ssm_smooth_eps`
- `p`: The length of the observation vectors, $\vec{y}_t$.


**returns** A $p \times p$ matrix with $\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$

Extract $\Var(\varepsilon_t|\vec{y}_{1:n})$ from vectors returned by `ssm_smooth_eps`


```stan

matrix ssm_smooth_eps_get_var(vector x, int p) {
  matrix[p, p] eps_var;
  eps_var = vector_to_symmat(x[(p + 1): ], p);
  return eps_var;
}


```



### ssm_smooth_eps

**Arguments**

- `filter`: Results of `ssm_filter`
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.


**returns** An array of vectors constaining $\hat{\vec{\varepsilon}}_t$ and $\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$ in the format described below.

The observation disturbance smoother

This calculates the mean and variance of the observation disturbances, $\vec{\varepsilon}_t$,
given the entire sequence, $\vec{y}_{1:n}$.




For Z`, `H`, T`, the array can have a size of 1, if it is not time-varying, or a size of $n$ (for `Z`, `H`) or $n - 1$ (for `T`),
if it is time varying.

The vectors returned by this function have $p + p (p + 1) / 2$ elements in this format,
$$
(\hat{\vec{\varepsilon}}_t', \VEC(\Var(\vec{\varepsilon}_t | \vec{y}_{1:n}))' )'
$$

value                                            length          start                 end
------------------------------------------------ --------------- --------------------- --------------------
$\hat{\vec{\varepsilon}}_t$                   $p$             $1$                 $p$
$\Var(\vec{\varepsilon}_t | \vec{y}_{1:n})$  $p (p + 1) / 2$ $p + 1$            $p + p (p + 1) / 2$


See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]


```stan

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
      L = ssm_update_L(Z_t, T_t, K);
      // r_{t - 1} and N_{t - 1}
      r = ssm_update_r(r, Z_t, v, Finv, L);
      N = ssm_update_N(N, Z_t, Finv, L);
      // eps_t and V(eps_t|y)
      eps = H_t * (Finv * v - K ' * r);
      var_eps = to_symmetric_matrix(H_t - H_t * (Finv + quad_form_sym(N, K)) * H_t);
      // saving
      res[t, :p] = eps;
      res[t, (p + 1): ] = symmat_to_vector(var_eps);
    }
  }
  return res;
}


```



### ssm_smooth_eta_size

**Arguments**

- `q`: The size of the state disurbance vector, $\vec{\eta}_t$.


**returns** The size of the vectors is $q + q (q + 1) / 2$.

The size of the vectors returned by `ssm_smooth_eta`


```stan

int ssm_smooth_eta_size(int q) {
  int sz;
  sz = q + symmat_size(q);
  return sz;
}


```



### ssm_smooth_eta_get_mean

**Arguments**

- `x`: A vector returned by `ssm_smooth_eta`
- `q`: The number of state disturbances, $\vec{\eta}_t$.


**returns** A $q \times 1$ vector with $\hat{\vec{\eta}}_t$.

Extract $\hat{\vec{\varepsilon}}_t$ from vectors returned by `ssm_smooth_eta`


```stan

vector ssm_smooth_eta_get_mean(vector x, int q) {
  vector[q] eta;
  eta = x[ :q];
  return eta;
}


```



### ssm_smooth_eta_get_var

**Arguments**

- `x`: A vector returned by `ssm_smooth_eta`
- `q`: The number of state disturbances, $\vec{\eta}_t$.


**returns** A $q \times q$ matrix with $\Var(\vec{\eta}_t | \vec{y}_{1:n})$.

Extract $\Var(\eta_t|\vec{y}_{1:n})$ from vectors returned by `ssm_smooth_eta`


```stan

matrix ssm_smooth_eta_get_var(vector x, int q) {
  matrix[q, q] eta_var;
  eta_var = vector_to_symmat(x[(q + 1): ], q);
  return eta_var;
}


```



### ssm_smooth_eta

**Arguments**

- `filter`: Results of `ssm_filter`
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.


**returns** An array of vectors constaining $\hat{\vec{\eta}}_t$ and $\Var(\vec{\eta}_t | \vec{y}_{1:n})$ in the format described below.

The state disturbance smoother

This calculates the mean and variance of the observation disturbances, $\vec{\eta}_t$,
given the entire sequence, $\vec{y}_{1:n}$.


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


```stan

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
      L = ssm_update_L(Z_t, T_t, K);
      r = ssm_update_r(r, Z_t, v, Finv, L);
      N = ssm_update_N(N, Z_t, Finv, L);
      eta = Q_t * R_t ' * r;
      var_eta = to_symmetric_matrix(Q_t - Q_t * quad_form_sym(N, R_t) * Q_t);
      // saving
      res[t, :q] = eta;
      res[t, (q + 1): ] = symmat_to_vector(var_eta);
    }
  }
  return res;
}



```



### ssm_smooth_state_mean

**Arguments**

- `filter`: The results of `ssm_filter`
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.


**returns** An array of size $n$ of $m \times 1$ vectors containing $\hat{\vec{\alpha}}_t$.

The fast state smoother

The fast state smoother calculates $\hat{\vec{\alpha}}_t = \E(\vec{\alpha}_t | \vec{y}_{1:n})$.
$$
\hat{\vec{\alpha}}_{t + 1} = \mat{T}_t \hat{\vec{\alpha}}_{t} + \mat{R}_t \mat{Q}_t \mat{R}'_t \vec{r}_t ,
$$
where $r_t$ is calcualted from the state disturbance smoother.
The smoother is initialized at $t = 1$ with $\hat{\vec{\alpha}}_t = \vec{a}_1 + \mat{P}_1 \vec{r}_0$.

Unlike the normal state smoother, it does not calculate the variances of the smoothed state.


For  `Z`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `Z`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

See [@DurbinKoopman2012, Sec 4.5.3 (eq 4.69)]


```stan

vector[] ssm_smooth_state_mean(vector[] filter,
                              matrix[] Z, vector[] c,
                              matrix[] T, matrix[] R, matrix[] Q) {
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
      L = ssm_update_L(Z_t, T_t, K);
      // r_{t - 1}
      r[t] = ssm_update_r(r[t + 1], Z_t, v, Finv, L);
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





```

##  Simulators and Smoothing Simulators




### ssm_sim_idx

**Arguments**

- `m`: The number of states
- `p`: The length of the observation vector
- `q`: The number of state disturbances


**returns** A 4 x 3 array of integers

Indexes of each component of `ssm_sim_rng` results.


The returned array has columns (length, start location, and end location)
for rows: $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\varepsilon}_t$, and $\vec{\eta}_t$ in the results of `ssm_sim_rng`.


```stan

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


```



### ssm_sim_size

**Arguments**

- `m`: The number of states
- `p`: The length of the observation vector
- `q`: The number of state disturbances


**returns** The number of elements

The number of elements in vectors returned by `ssm_sim_rng` results.


```stan

int ssm_sim_size(int m, int p, int q) {
  int sz;
  sz = ssm_sim_idx(m, p, q)[4, 3];
  return sz;
}


```



### ssm_sim_get_y

**Arguments**

- `x`: vector of results from `ssm_sim_rng`
- `m`: The number of states
- `p`: The length of the observation vector
- `q`: The number of state disturbances


**returns** vector A $p \times 1$ vector with $\vec{y}_t$.

Extract $\vec{y}_t$ from vectors returned by `ssm_sim_rng`.


```stan

vector ssm_sim_get_y(vector x, int m, int p, int q) {
  vector[p] y;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  y = x[idx[1, 2]:idx[1, 3]];
  return y;
}


```



### ssm_sim_get_a

**Arguments**

- `x`: vector of results from `ssm_sim_rng`
- `m`: The number of states
- `p`: The length of the observation vector
- `q`: The number of state disturbances


**returns** A $m \times 1$ vector with $\vec{\alpha}_t$.

Extract $\vec{\alpha}_t$ from vectors returne by `ssm_sim_rng`.


```stan

vector ssm_sim_get_a(vector x, int m, int p, int q) {
  vector[m] a;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  a = x[idx[2, 2]:idx[2, 3]];
  return a;
}


```



### ssm_sim_get_eps

**Arguments**

- `x`: vector of results from `ssm_sim_rng`
- `m`: The number of states
- `p`: The length of the observation vector
- `q`: The number of state disturbances


**returns** vector A $p \times 1$ vector with $\vec{\varepsilon}_t$.

Extract $\vec{\varepsilon}_t$ from vectors returne by `ssm_sim_rng`.


```stan

vector ssm_sim_get_eps(vector x, int m, int p, int q) {
  vector[p] eps;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  eps = x[idx[3, 2]:idx[3, 3]];
  return eps;
}


```



### ssm_sim_get_eta

**Arguments**

- `x`: vector of results from `ssm_sim_rng`
- `m`: The number of states
- `p`: The length of the observation vector
- `q`: The number of state disturbances


**returns** vector A $q \times 1$ vector with $\vec{\eta}_t$.

Extract $\vec{\eta}_t$ from vectors returne by `ssm_sim_rng`.


```stan

vector ssm_sim_get_eta(vector x, int m, int p, int q) {
  vector[q] eta;
  int idx[4, 3];
  idx = ssm_sim_idx(m, p, q);
  eta = x[idx[4, 2]:idx[4, 3]];
  return eta;
}


```



### ssm_sim_rng

**Arguments**

- `n`: Number of time observations to draw, $t = 1, \dots, n$.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** of size $n$ of vectors with Draw $\vec{y}_t$, $\vec{\alpha}_t$, $\vec{\eta}_t$ and $\vec{\varepsilon}_t$. See the description.

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


```stan

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
    int idx[4, 3];

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
      ret[t, idx[3, 2]:idx[3, 3]] = eps;
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
      } else {
        // don't forecast alpha_{t + 1}, so don't draw eta_t
        eta = zero_q;
      }
      // save eta_t; alpha is saved at the start of the loop.
      ret[t, idx[4, 2]:idx[4, 3]] = eta;
    }
  }
  return ret;
}


```


##  Simulation Smoothers




### ssm_simsmo_states_rng

**Arguments**

- `filter`: A length $n$ array with results from `ssm_filter`.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** Array of size $n$ of $m \times 1$ vectors containing a single draw from $(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.

State simulation smoother

Draw samples from the posterior distribution of the states,
$\tilde{\vec{\alpha}}_{1:n} \sim p(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.


For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].


```stan

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
      alpha_hat = ssm_smooth_state_mean(filter, Z, c, T, R, Q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter with simulated y's
      filter_plus = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      alpha_hat_plus = ssm_smooth_state_mean(filter_plus, Z, c, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha_hat[i]);
      }
    }
    return draws;
}



```



### ssm_simsmo_states_miss_rng

**Arguments**

- `filter`: A length $n$ array with results from `ssm_filter`.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
- `p_t`: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- `y_idx`: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** Array of size $n$ of $m \times 1$ vectors containing a single draw from $(\vec{\alpha}_{1:n} | \vec{y}_{1:n})$.

State simulation smoother, as in `ssm_simsmo_states_rng`, allowing for missing values.


```stan

vector[] ssm_simsmo_states_miss_rng(vector[] filter,
                      vector[] d, matrix[] Z, matrix[] H,
                      vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                      vector a1, matrix P1, int[] p_t, int[,] y_idx) {
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
      alpha_hat = ssm_smooth_state_mean(filter, Z, c, T, R, Q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter with simulated y's
      filter_plus = ssm_filter_miss(y, d, Z, H, c, T, R, Q, a1, P1, p_t, y_idx);
      // mean correct epsilon samples
      alpha_hat_plus = ssm_smooth_state_mean(filter_plus, Z, c, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha_hat[i]);
      }
    }
    return draws;
}


```



### ssm_simsmo_eta_rng

**Arguments**

- `filter`: A length $n$ array with results from `ssm_filter`.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** Array of size $n$ of $q \times 1$ vectors containing a single draw from $(\vec{\eta}_{1:n} | \vec{y}_{1:n})$.

State disturbance simulation smoother

Draw samples from the posterior distribution of the observation disturbances,
$\tilde{\vec{\eta}}_{1:n} \sim p(\vec{\eta}_{1:n} | \vec{y}_{1:n})$.


For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].


```stan

vector[] ssm_simsmo_eta_rng(vector[] filter,
                            vector[] d, matrix[] Z, matrix[] H,
                            vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                            vector a1, matrix P1) {
    vector[dims(Q)[2]] draws[size(filter)];
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
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eta_size(q)] eta_hat[n];
      vector[ssm_smooth_eta_size(q)] eta_hat_plus[n];
      // get smoothed etas
      eta_hat = ssm_smooth_eta(filter, Z, T, R, Q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter_plus = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct eta samples
      eta_hat_plus = ssm_smooth_eta(filter_plus, Z, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_eta(sims[i], m, p, q)
                    - ssm_smooth_eta_get_mean(eta_hat_plus[i], q)
                    + ssm_smooth_eta_get_mean(eta_hat[i], q));
      }
    }
    return draws;
}



```



### ssm_simsmo_eta_miss_rng

**Arguments**

- `filter`: A length $n$ array with results from `ssm_filter`.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
- `p_t`: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- `y_idx`: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** Array of size $n$ of $q \times 1$ vectors containing a single draw from $(\vec{\eta}_{1:n} | \vec{y}_{1:n})$.

State disturbance simulation smoother, as in `ssm_simsmo_eta_rng`, but allowing for missing values in $\vec{y}_t$.


```stan

vector[] ssm_simsmo_eta_miss_rng(vector[] filter,
                            vector[] d, matrix[] Z, matrix[] H,
                            vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                            vector a1, matrix P1, int[] p_t, int[,] y_idx) {
    vector[dims(Q)[2]] draws[size(filter)];
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
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eta_size(q)] eta_hat[n];
      vector[ssm_smooth_eta_size(q)] eta_hat_plus[n];
      // get smoothed etas
      eta_hat = ssm_smooth_eta(filter, Z, T, R, Q);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter_plus = ssm_filter_miss(y, d, Z, H, c, T, R, Q, a1, P1, p_t, y_idx);
      // mean correct eta samples
      eta_hat_plus = ssm_smooth_eta(filter_plus, Z, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_eta(sims[i], m, p, q)
                    - ssm_smooth_eta_get_mean(eta_hat_plus[i], q)
                    + ssm_smooth_eta_get_mean(eta_hat[i], q));
      }
    }
    return draws;
}


```



### ssm_simsmo_eps_rng

**Arguments**

- `filter`: A length $n$ array with results from `ssm_filter`.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.


**returns** Array of size $n$ of $p \times 1$ vectors containing a single draw from $(\vec{\varepsilon}_{1:n} | \vec{y}_{1:n})$.

Observation disturbance simulation smoother.

Draw samples from the posterior distribution of the observation disturbances,
$\tilde{\vec{\varepsilon}}_{1:n} \sim p(\vec{\varepsilon}_{1:n} | \vec{y}_{1:n})$.


For `d`, `Z`, `H`, `c`, `T`, `R`, `Q` the array can have a size of 1, if it is
not time-varying, or a size of $n$ (for `d`, `Z`, `H`) or $n - 1$ (for `c`, `T`, `R`, `Q`)
if it is time varying.

This draws samples using mean-correction simulation smoother of [@DurbinKoopman2002].
See [@DurbinKoopman2012, Sec 4.9].


```stan

vector[] ssm_simsmo_eps_rng(vector[] filter,
                      vector[] d, matrix[] Z, matrix[] H,
                      vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                      vector a1, matrix P1) {
    vector[dims(Z)[2]] draws[size(filter)];
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
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eps_size(p)] eps_hat_plus[n];
      vector[ssm_smooth_eps_size(p)] eps_hat[n];

      // get smoothed values of epsilon
      eps_hat = ssm_smooth_eps(filter, Z, H, T);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter_plus = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      // mean correct epsilon samples
      eps_hat_plus = ssm_smooth_eps(filter_plus, Z, H, T);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_eps(sims[i], m, p, q)
                    - ssm_smooth_eps_get_mean(eps_hat_plus[i], p)
                    + ssm_smooth_eps_get_mean(eps_hat[i], p));
      }
    }
    return draws;
}



```



### ssm_simsmo_eps_miss_rng

**Arguments**

- `filter`: A length $n$ array with results from `ssm_filter`.
- `d`: Observation intercept, $\vec{d}_t$. An array of $p \times 1$ vectors.
- `Z`: Design matrix, $\mat{Z}_t$. An array of $p \times m$ matrices.
- `H`: Observation covariance matrix, $\mat{H}_t$. An array of $p \times p$ matrices.
- `c`: State intercept, $\vec{c}_t$. An array of $m \times 1$ vectors.
- `T`: Transition matrix, $\mat{T}_t$. An array of $m \times m$ matrices.
- `R`: State covariance selection matrix, $\mat{R} _t$. An array of $p \times q$ matrices.
- `Q`: State covariance matrix, $\mat{Q}_t$. An array of $q \times q$ matrices.
- `a1`: Expected value of the intial state, $a_1 = \E(\alpha_1)$. An $m \times 1$ matrix.
- `P1`: Variance of the initial state, $P_1 = \Var(\alpha_1)$. An $m \times m$ matrix.
- `p_t`: A length $n$ array with the number of non-missing elements in the observation vector, $\vec{y}_t$, at each $t \in 1, \dots, n$.
- `y_idx`: A length $n \times p$ array of integers. The first $p_t$ elments of this array indexes of thenon-missing values of $y$. Elements $1:p_t$ should be between $1$ and $p$; elements $p_t:p$ are zero, and are not used.


**returns** Array of size $n$ of $p \times 1$ vectors containing a single draw from $(\vec{\varepsilon}_{1:n} | \vec{y}_{1:n})$.

Observation disturbance simulation smoother, as in `ssm_simsmo_eps_rng`, but allowing for missing values in $\vec{y}_t$.


```stan
vector[] ssm_simsmo_eps_miss_rng(vector[] filter,
                      vector[] d, matrix[] Z, matrix[] H,
                      vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                      vector a1, matrix P1, int[] p_t, int[,] y_idx) {
    vector[dims(Z)[2]] draws[size(filter)];
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
      vector[p] y[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[ssm_smooth_eps_size(p)] eps_hat_plus[n];
      vector[ssm_smooth_eps_size(p)] eps_hat[n];

      // get smoothed values of epsilon
      eps_hat = ssm_smooth_eps(filter, Z, H, T);
      // simulate unconditional disturbances and observations
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      // filter simulated y's
      filter_plus = ssm_filter_miss(y, d, Z, H, c, T, R, Q, a1, P1, p_t, y_idx);
      // mean correct epsilon samples
      eps_hat_plus = ssm_smooth_eps(filter_plus, Z, H, T);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_eps(sims[i], m, p, q)
                    - ssm_smooth_eps_get_mean(eps_hat_plus[i], p)
                    + ssm_smooth_eps_get_mean(eps_hat[i], p));
      }
    }
    return draws;
}


```


##  Stationary




### pacf_to_acf

**Arguments**

- `x`: A vector of coefficients of a partial autocorrelation function


**returns** A vector of coefficients of an Autocorrelation function

Partial Autocorrelations to Autocorrelations


```stan

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


```



### constrain_stationary

**Arguments**

- `x`: An unconstrained vector in $(-\infty, \infty)$


**returns** A vector of coefficients for a stationary AR or inverible MA process.

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


```stan

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




```



### acf_to_pacf

**Arguments**

- `x`: Coeffcients of an autocorrelation function.


**returns** A vector of coefficients of the corresponding partial autocorrelation function.

Convert coefficients of an autocorrelation function to partial autocorrelations.


```stan

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


```



### unconstrain_stationary

**Arguments**

- `x`: Coeffcients of an autocorrelation function.


**returns** Coefficients of the corresponding partial autocorrelation function.

Transform from stationary and invertible space to $(-\infty, \infty)$.


```stan

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




```



### kronecker_prod

**Arguments**

- `A`: An $m \times n$ matrix
- `B`: A $p \times q$ matrix


**returns** An $mp \times nq$ matrix.

Kronecker product

The Kronecker product of a $\mat{A}$ and $\mat{B}$ is
$$
\mat{A} \otimes \mat{B} =
\begin{bmatrix}
a_{11} \mat{B} \cdots a_{1n} \mat{B} \
\vdots & \ddots & vdots \
a_{m1} \mat{B} & \cdots & a_{mn} \mat{B}
\end{bmatrix} .
$$


```stan

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


```



### stationary_cov

**Arguments**

- `T`: The $m \times m$ transition matrix
- `RQR`: The $m \times m$ system covarariance matrix, $\mat{R} \mat{Q} \mat{R}\T$.


**returns** An $m \times m$ matrix with the stationary covariance matrix.

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


```stan

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

```
