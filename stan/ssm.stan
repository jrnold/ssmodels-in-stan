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

// // Partial Autocorrelations to Autocorrelations transformation
// // Maps (-1, 1)^p to AR space
// // Translated from R stats C function partrans
// // https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/pacf.c
// vector pacf_to_acf(vector x) {
//   vector[size(x)] ret;
//   vector[size(x)] work;
//   real a;
//   int p;
//   p = size(x);
//
//   // elements must be between -1 and 1
//   for (i in 1:p) {
//     if (x[i] > 1) {
//       reject("x is greater than 1");
//     } else if (x[i] < -1) {
//       reject("x is less than -1");
//     }
//   }
//
//   work = x;
//   ret = work;
//   /* run the Durbin-Levinson recursions to find phi_{j.},
//      j = 2, ..., p and phi_{p.} are the autoregression coefficients */
//   for(j = 2:p) {
//   	a = ret[j];
//     for(k in 1:j) {
//       // TODO: possible off by one error here
//       work[k] = work[k] - a * ret[j - k - 1];
//     }
//     ret = work;
//   }
//   return ret
// }
//
// // Autocorrelations to Partial Autocorrelations
// // Maps AR space to (-1, 1)^p
// // Translated from R stats C function invpartrans
// // https://github.com/wch/r-source/blob/e5b21d0397c607883ff25cca379687b86933d730/src/library/stats/src/pacf.c
// vector acf_to_pacf(vector x) {
//     real a;
//     vector[size(x)] ret;
//     vector[size(x)] work;
//     int p;
//
//     p = size(x);
//     work = x;
//     for(j in 1:p) {
//       work[j] = phi[j];
//       ret[j] = work[j];
//     }
//     /* Run the Durbin-Levinson recursions backwards
//        to find the PACF phi_{j.} from the autoregression coefficients */
//     for(i in 2:p) {
//       int j;
//       j = p - i;
// 	    a = ret[j];
// 	    for(k in 1:j) {
//         // TODO: possible off by one error here
//         work[k]  = (ret[k] + a * ret[j - k - 1]) / (1 - a * a);
//       }
//       ret = work;
//     }
//     // let user check that outputs are between -1 and 1.
//     return ret;
// }
