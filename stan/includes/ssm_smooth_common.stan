///
/// # Common Smoother Functions
///

/**

$$
\vec{r}_t = \mat{Z}' \mat{F}^{-1}_t \vec{v}_t + \mat{L}' \vec{r}_{t - 1}
$$

See [@DurbinKoopman2012, p. 91]
*/
vector ssm_smooth_update_r(vector r, matrix Z, vector v, matrix Finv,
                           matrix L) {
  vector[num_elements(r)] r_new;
  r_new = Z ' * Finv * v + L ' * r;
  return r_new;
}

/**

$$
\mat{N}_t = \mat{Z})_t' \mat{F}^{-1}_t \mat{Z}_t + \mat{L}_t' \mat{N}_t \mat{L}_t
$$

See [@DurbinKoopman2012, p. 91]
*/
matrix ssm_smooth_update_N(matrix N, matrix Z, matrix Finv, matrix L) {
  matrix[rows(N), cols(N)] N_new;
  N_new = quad_form(Finv, Z) + quad_form(N, L);
  return N_new;
}
