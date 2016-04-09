// -*- mode: stan -*-
int int_number_true(int x) {
  int n;
  n <- 0;
  for (i in 1:size(x)) {
    if (int_step(x[i])) {
       n <- n + 1;
    }
  }
  return n;
}

// symmetric matrix to vector
vector symmat_to_vector(matrix X) {
  vector[(rows(X) * (rows(X) - 1)) / 2] v;
  int k;
  k <- 1;
  for (i in 1:rows(X)) {
    for (j in 1:i) {
      v[k] <- X[i, j];
      k <- k + 1;
    }
  }
  return v;
}

// matrix to symmetric vector
matrix vector_to_symmat(vector v) {
  matrix[floor(sqrt(size(v))), floor(sqrt(size(v)))] X;
  k <- 1;
  for (i in 1:rows(X)) {
    for (j in 1:i) {
      if (i != j) {
        X[i, j] <- X[j, i] <- v[k];
      } else {
        X[i, i] <- v[k];
      }
      k <- k + 1;
    }
  }
  return X;
}

matrix select_vec(int x, int i) {
  vector[size(i)] y;
  for (i in 1:size(i)) {
    y[i] <- x[i]
  }
  return y;
}
matrix select_mat_rc(int X, int indices) {
  matrix[size(indices), size(indices)] Y;
  for (j in 1:size(indices)) {
    for (k in 1:size(indices)) {
      Y[j, k] <- X[i[j], i[k]];
    }
  }
  return Y;
}
matrix select_mat_row(int X, int indices) {
  matrix[size(indices), cols(X)] Y;
  for (i in 1:size(indices)) {
    for (j in 1:cols(X)) {
      Y[i, j] <- X[indices[i], j];
    }
  }
  return Y;
}

matrix select_mat_col(int X, int indices) {
  matrix[rows(X), size(indices)] Y;
  for (i in 1:rows(X)) {
    for (j in 1:size(indices)) {
      Y[i, j] <- X[i, indices[j]];
    }
  }
  return Y;
}

// select with masks
vector select_vec_mask(vector x, int mask) {
  vector[sum(mask)] y;
  int k;
  k <- 1;
  for (i in 1:size(mask)) {
    if (int_step(mask)) {
      y[k] <- x[i];
      k <- k + 1;
    }
  }
  return y;
}

matrix make_symmteric_matrix(matrix X) {
  return 0.5 * (X + X');
}

// kalman_filter_predict_obs
// kalman_filter_predict_state
// kalman_filter_update_state
// kalman_filter_log
// kalman_filter_backwards_sample
}
