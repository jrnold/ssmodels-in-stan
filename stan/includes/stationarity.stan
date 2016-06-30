/**
Partial Autocorrelations to Autocorrelations

@param vector of coefficients of a partial autocorrelation function
@return vector of coefficients of an Autocorrelation function

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

/**
Constrain vector of coefficients to the stationary and intertible region for AR or MA functions.

@param vector x An unconstrained vector in (-Inf, Inf)
@return vector y A vector of coefficients for a stationary AR or inverible MA process.

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

/**
Convert coefficients of an autocorrelation function to partial autocorrelations.

@param vector x Coeffcients of an autocorrelation function.
@return vector y Coefficients of the corresponding partial autocorrelation function.

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

/**
Transform from stationary and invertible space to (-Inf, Inf)

@param vector x Coeffcients of an autocorrelation function.
@return vector y Coefficients of the corresponding partial autocorrelation function.

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

/**
Kronecker product

The Kronecker product of a $A$ and $B$ is
$$
A \crossprod B =
\begin{bmatrix}
a_{11} B \cdots a_{1n} B \\
\vdots & \ddots & vdots \\
a_{m1} B & \cdots & a_{mn} B
\end{bmatrix} .
$$

@param matrix A An $m \times n$ matrix
@param matrix B A $p \times q$ matrix
@return A $mp \times nq$ matrix.

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
      col_end = (j - 1) * q + 1;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}
re

/**
Find the covariance of the stationary distribution of an ARMA model

@param matrix T The $m \times m$ transition matrix
@param matrix R The $m \times q$ system disturbance selection matrix
@param A $m \times m$ matrix with the stationary covariance matrix.

The initial conditions are $\alpha_1 \sim N(0, \sigma^2 Q_0),
where $Q_0$ is the solution to
$$
(T \crossprod T) vec(Q_0) = vec(R R')
$$
where $vec(Q_0)$ and $vec(R R')$ are the stacked columns of $Q_0$ and $R R'$

See [@DurbinKoopmans2012, Sec 5.6.2].

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
