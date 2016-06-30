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
  i = 0;
  while (n > 0) {
    i = i + 1;
    n = n - 1;
  }
  return i;
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
