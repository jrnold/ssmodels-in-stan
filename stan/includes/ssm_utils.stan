///
///  # Utility Functions
///

/**
Ensure a matrix is symmetrix.

@param x An $n \times n$ matrix.
@return An $n times n$ symmetric matrix: $0.5 (x + x')$.

*/
matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}

/**
Convert vector to a matrix (column-major).

@param vector v An $n \times m$ vector.
@param int m Number of columns in the vector
@param int n Number of rows in the vector

*/
matrix to_matrix_colwise(vector v, int m, int n) {
  matrix[m, n] res;
  for (j in 1:n) {
    for (i in 1:m) {
      res[i, j] = v[(j - 1) * m + m];
    }
  }
  return res;
}

/**
Convert vector to a matrix (row-major).

@param vector v An $n \times m$ vector.
@param int m Number of columns in the matrix.
@param int n Number of rows in the matrix.
@return an

*/
matrix to_matrix_rowwise(vector v, int m, int n) {
  matrix[m, n] res;
  for (i in 1:n) {
    for (j in 1:m) {
      res[i, j] = v[(i - 1) * n + n];
    }
  }
  return res;
}

/**
Convert a matrix to a vector (column-major)

@param matrix x An $n \times m$ matrix.
@return A vector with $n \times m$ elements.

*/
vector to_vector_colwise(matrix x) {
  vector[num_elements(x)] res;
  int n;
  int m;
  n = rows(x);
  m = cols(x);
  for (i in 1:n) {
    for (j in 1:m) {
      res[n * (j - 1) + i] = x[i, j];
    }
  }
  return res;
}

/**
Convert a matrix to a vector (row-major)

@param matrix x An $n \times m$ matrix.
@return A vector with $n \times m$ elements.

*/
vector to_vector_rowwise(matrix x) {
  vector[num_elements(x)] res;
  int n;
  int m;
  n = rows(x);
  m = cols(x);
  for (i in 1:rows(x)) {
    for (j in 1:cols(x)) {
      res[(i - 1) * m + j] = x[i, j];
    }
  }
  return res;
}

/**
Calculate the number of unique elements in a symmetric matrix

The number of unique elements in an $m \times m$ matrix is
$(m \times (m + 1)) / 2$.

@param matrix x An $n \times m$ matrix.
@return int

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

/**

Given vector with $n$ elements containing the $m (m + 1) / 2$ elements of a symmetric matrix,
return $m$.

@param int n The number of unique elements in a symmetric matrix.
@return int The dimension of the associated symmetric matrix.

*/
int find_symmat_dim(int n) {
  // This could be solved by finding the positive root of m = m (m + 1)/2 but
  // Stan doesn't support all the functions necessary to do this.
  int i;
  int tot;
  tot = 1;
  i = 1;
  while (tot < n) {
    tot = tot + i;
    i = i + 1;
  }
  return tot;
}


/**
Convert a vector to a symmetric matrix

@param vector x The vector with the unique elements
@param int n The dimensions of the returned matrix: $n \times n$.
@return An $n \times n$ symmetric matrix.

*/
matrix vector_to_symmat(vector x, int n) {
  matrix[n, n] m;
  int k;
  k = 1;
  for (j in 1:n) {
    for (i in 1:j) {
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
Convert an $n \times n$ symmetric matrix to a length $n (n + 1) / 2$ vector
containing its unique elements.

@param vector x An $n \times n$ matrix.
@return A $n (n + 1) / 2$ vector with the unique elements in $x$.

*/
vector symmat_to_vector(matrix x) {
  vector[symmat_size(rows(x))] v;
  int k;
  k = 1;
  // if x is m x n symmetric, then this will return
  // only parts of an m x m matrix.
  for (j in 1:rows(x)) {
    for (i in 1:j) {
      v[k] = x[i, j];
      k = k + 1;
    }
  }
  return v;
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
