
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
