// -*- mode: stan -*-

// Utility Functions

matrix make_symmetric_matrix(matrix x) {
  return 0.5 * (x + x');
}

matrix to_matrix_colwise(vector v, int r, int c) {
  matrix[r, c] res;
  for (i in 1:size(v)) {
      res[(i - 1) % r + 1, (i - 1) / r + 1] <- v[i]
  }
  return res;
}

matrix to_matrix_rowwise(vector v, int r, int c) {
  matrix[r, c] res;
  for (i in 0:(size(v) - 1)) {
      res[(i - 1) / c + 1, (i - 1) % c + 1] <- v[i]
  }
  return res;
}

matrix to_vector_colwise(matrix x) {
  vector[num_elements(x)] res;
  for (r in 1:rows(x)) {
    for (c in 1:cols(x)) {
      res[r * (c - 1) + r] <- x[r, c]
    }
  }
  return res;
}

matrix to_vector_rowwise(matrix x) {
  vector[num_elements(x)] res;
  for (r in 1:rows(x)) {
    for (c in 1:cols(x)) {
      res[(r - 1) * c + c] <- x[r, c]
    }
  }
  return res;
}


// Filtering

vector ssm_filter_pred_obs(vector y, matrix Z, matrix H, vector a, matrix P) {

}
